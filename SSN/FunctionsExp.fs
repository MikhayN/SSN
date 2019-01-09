module FunctionsExp

open Types

open FSharp.Stats
open FSharpAux

open FSharp.Plotly

open System

open PQ

module SSN =
    ////// Prepare data

    /// get id list from protein list
    let groupIDFn (items : Item []) = 
        items |> Array.map (fun ii -> ii.ID) 

    ////// Calculation

    /// calculate euclidian distance between two vectors with weight (optional)
    let weightedEuclidean (weightL: seq<float> option) v1 v2 = 
            let n = 
                v1
                |> Seq.length
            let weightL' =
                match weightL with
                |Some x -> x
                |None -> Seq.initRepeatValue n 1.
            Seq.zip3 weightL' v1 v2
            |> Seq.fold (fun d (w12,e1,e2) -> d + w12*((e1 - e2) * (e1 - e2))) 0.
            |> sqrt


    /// apply weighted Euclidean distance to create a matrix of two vectors of data
    let distanceMatrixWeighted weightL (data:float[][]) =
        let m = Array2D.zeroCreate (data.Length) (data.Length)
        for rowI in 0..data.Length-1 do
            for colI in 0..rowI do
                let tmp = weightedEuclidean weightL data.[rowI] data.[colI] 
                m.[colI,rowI] <- tmp
                m.[rowI,colI] <- tmp
        m

    /// calculate matrix for a list of Item
    let distMatrixWeightedOf f weightL (itemOutList: array<Types.Item>) =
        let listData = 
            itemOutList
            |> Array.map (fun i -> i.dataL) 
        if (listData.Length > 1) then 
                f weightL listData 
        else 
                Array2D.zeroCreate 2 2

    // change here if want smth not MAX
    /// find the max of dissimilarities for a protein (idCurrent) in a group of Item (idGroup) by their ids


    // Calculation of gain, given protein kinetics and kinetics of all Item in current and predeccesor node 

    /// search through matrix and find max dissimilarity for each element in itemsToSum and sum them
    let dSumFn itemsToSum itemsWhereToFind matrix =
        
        let findMaxDist idCurrent (idGroup: int array) (matrix: float [,]) = 
            idGroup
            |> Array.fold (fun maxSoFar i -> max maxSoFar matrix.[i,idCurrent]) 0.
        
        itemsToSum
        |> Array.fold (fun acc i -> acc + (findMaxDist i itemsWhereToFind matrix)) 0.

    /// calculate step gain between parent and child nodes given the root number with given function for G_s
    let getStepGainFn fn itemsChild itemsParent (nRoot: int) matrix : float =

        let dSumF itemsToSum itemsWhereToFind matrix =
        
            let findMaxDistInMatrix idCurrent (idGroup: int array) (matrix: float [,]) = 
                idGroup
                |> Array.fold (fun maxSoFar i -> max maxSoFar matrix.[i,idCurrent]) 0.
        
            itemsToSum
            |> Array.fold (fun acc i -> acc + (findMaxDistInMatrix i itemsWhereToFind matrix)) 0.

        let dPredSum = dSumF itemsChild itemsParent matrix
        let dCurrSum = dSumF itemsChild itemsChild matrix
        fn dCurrSum dPredSum itemsChild.Length nRoot 

    let getStepGainFnPQ fn itemsChild itemsParent (nRoot: int) (pq: MaxIndexPriorityQueue<float> []) : float =
        
        itemsChild |> Array.iter (fun i -> pq.[i].LeaveGroup itemsParent)
        let dPredSum = itemsChild |> Array.sumBy (fun x -> pq.[x].Top())

        itemsChild |> Array.iter (fun i -> pq.[i].LeaveGroup itemsChild)
        let dCurrSum = itemsChild |> Array.sumBy (fun x -> pq.[x].Top())

        fn dCurrSum dPredSum itemsChild.Length nRoot 

    /// function to calculate configuration gain
    let confGainFn (items: Map<string,Types.Node<string,Types.Item>>) =                              /// complex input
        items
        |> Map.fold (fun state key value -> state+(value.GroupGain)) 0. 

    ////// Tree structure

    /// break group in terms find all SO-groups of items (not a single one!) in a current node, 
    /// looking at path list and group by item.next_depth and also find Item, that stays left in current node
    let rec breakGroup (items: Types.Item array) depth =
    
        if (items |> Array.forall (fun ii -> ii.OriginalBin.Length>depth+1))
            && (items |> Array.forall (fun ii -> ii.OriginalBin.[depth+1]=items.[0].OriginalBin.[depth+1])) then // checking for the tunnels
            let newItems = 
                items
                |> Array.map (fun x -> {x with OriginalBin=(x.OriginalBin |> Array.removeIndex (depth))})
            breakGroup newItems depth
        else
            items 
            |> Array.map (fun i -> 
                if i.OriginalBin.Length<=(depth+1) then 
                     {i with 
                        OriginalBin = Array.append i.OriginalBin [|sprintf "%s%i" "p" i.ID|]; 
                        BinL = Array.append i.BinL [|sprintf "%s%i" "p" i.ID|]}
                else i)
            |> Array.groupBy (fun i -> i.OriginalBin.[depth+1])
            |> Map.ofArray

    /// Intermediate partition function: integer partition of n to k parts
    let schemeGenerator n k =
        let schemeArray = Array.init k (fun i -> 0)
        let rec loop nRest kRest prevSum maxValue =
            seq [
                if (prevSum=n) || (kRest=0) then
                    let temp = schemeArray.[0 .. (k-1)]             // to be able to write inside array
                    yield temp
                else
                    let half = int (ceil((float nRest)/(float kRest)))
                    let list = List.init ( (min maxValue (nRest-kRest+1)) - half + 1 ) (fun i -> half + i )
                    for a in list do 
                        schemeArray.[(k-kRest)] <- a
                        yield! loop (nRest-a) (kRest-1) (prevSum+a) a 
            ]
        loop n k 0 n

    ////// Create tree

    //
    ////// Add hierarchical clustering
    //



    type Mode =                                 // creating/alterating tree approach
        |MM         // MapMan with broken leaves
        |SSN        // with simplification overall

    /// create tree function with two modes: MM - just read and show original MapMan ontology;
    /// SSN - process the MM tree into optimal SSN structure
    /// gainFn - gain formula, 
    /// kmeanswapFn - function for swapping,
    /// clusterFn - function for pre-clustering in case of more than 50 singletons as leaves
    let createTree gainFn kmeanswapFn clusterFn (weight: seq<float> option) (mode: Mode) (rootGroup: Types.Item array) = 

        let nRoot = rootGroup.Length

        let matrix = 
            rootGroup
            |> distMatrixWeightedOf distanceMatrixWeighted weight

        // calculation for one node    
        let rec loop (nodeMembers: Types.Item array) depth dPredSum =
        
            /// sum of max dist within a node 
            let dCurrSum = dSumFn (groupIDFn nodeMembers) (groupIDFn nodeMembers) matrix

            /// to calc step from parent to current node
            let stepGain = (gainFn dCurrSum dPredSum nodeMembers.Length nRoot)

            let children = 
                match mode with
                |MM -> // MapMan with broken leaves
                    if nodeMembers.Length=1 then
                        Map.empty
                    else 
                        (breakGroup nodeMembers depth)
                        |> Map.map (fun key nodes -> 
                                        let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrix

                                        (loop nodes (depth+1) dPredSum'))

                |SSN -> // with simplification overall
                    if (nodeMembers.Length=1) 
                            || (nodeMembers.Length=0) 
                            || (String.contains "|" (String.Concat nodeMembers.[0].BinL))  then 
                        Map.empty
                    else 
                        (breakGroup nodeMembers depth) // get a map with original grouping as (subbin, Item)
                        |> (fun x -> 
                                    let singletons = x |> Map.filter (fun k pList -> pList.Length=1)
                                    if (singletons.Count>50) then
                                          let k = 50
                                          let list = singletons |> Map.toArray |> Array.map (snd) |> Array.concat |> List.ofArray
                                          let newClusters = clusterFn k weight list //mapOfSOM (snd (applySOM weight list 10 k)) list
                                          let oldClusters = x |> Map.filter (fun k pList -> pList.Length<>1)
                                          Map.fold (fun s k v -> Map.add k v s) newClusters oldClusters

                                     else 
                                          x)
                        |> kmeanswapFn gainFn matrix depth
                        |> Seq.fold (fun (singles,best) i -> 
                                            let newNodes = 
                                                if singles=Map.empty then
                                                    i
                                                    |> Map.fold (fun state key nodes ->
                                                                        let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrix
                                                                        state 
                                                                        |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)
                                                else
                                                    i
                                                    |> Map.fold (fun state key nodes ->
                                                                        match (singles.TryFind key) with
                                                                        | None ->
                                                                            let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrix
                                                                            state 
                                                                            |> Map.add key (loop nodes (depth+1) dPredSum')
                                                                        | Some x ->
                                                                            state 
                                                                            |> Map.add key x                                                                   
                                                                    ) (Map.empty)
                                                                    
                                            let best' =
                                                if (confGainFn newNodes) > (confGainFn best) then  // compare configuration gains to get the best
                                                    newNodes  
                                                else 
                                                    best
                                            if (singles = Map.empty) then
                                                (newNodes, best')
                                            else
                                                (singles, best')
                                            
                                        ) (Map.empty,Map.empty) // here as state should be this singles (first) saved and optimal conf
                        |> snd

            
            {
            StepGain = stepGain; 
            ConfGain = (confGainFn children, children |> Map.toList |> List.map fst);
            GroupGain = max stepGain (confGainFn children);
            Member = nodeMembers;
            Children = 
                if  (mode=MM) || ((confGainFn children)>stepGain) then 
                    children;
                else 
                    Map.empty
            }
    
        loop rootGroup 0 0.

    let createTreeReducing pre_k gainFn kmeanswapFn clusterFn (weight: seq<float> option) (mode: Mode) (rootGroup: Types.Item array) = 

        let nRoot = rootGroup.Length
        //let pre_k = 19//Array.create (rootGroup.[0].dataL.Length) 2 |> Array.reduce (fun a b -> a*b)

        let matrix = 
            rootGroup
            |> distMatrixWeightedOf distanceMatrixWeighted weight

        // calculation for one node    
        let rec loop (nodeMembers: Types.Item array) depth dPredSum =
        
            /// sum of max dist within a node 
            let dCurrSum = dSumFn (groupIDFn nodeMembers) (groupIDFn nodeMembers) matrix

            /// to calc step from parent to current node
            let stepGain = (gainFn dCurrSum dPredSum nodeMembers.Length nRoot)

            let children = 
                match mode with
                |MM -> // MapMan with broken leaves
                    if nodeMembers.Length=1 then
                        Map.empty
                    else 
                        (breakGroup nodeMembers depth)
                        |> Map.map (fun key nodes -> 
                                        let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrix

                                        (loop nodes (depth+1) dPredSum'))

                |SSN -> // with simplification overall
                    if (nodeMembers.Length=1) 
                            || (nodeMembers.Length=0) 
                            || (String.contains "|" (String.Concat nodeMembers.[0].BinL))  then 
                        Map.empty
                    else 
                        (breakGroup nodeMembers depth) // get a map with original grouping as (subbin, Item)
                        |> (fun x ->
                                    if (x.Count>pre_k) then
                                          let k = pre_k
                                          let newClusters = clusterFn k weight x //mapOfSOM (snd (applySOM weight list 10 k)) list
                                          let preG = newClusters |> Map.toArray |> Array.sumBy (fun (s,items) -> getStepGainFn gainFn (items |> groupIDFn) (nodeMembers |> groupIDFn) nodeMembers.Length matrix)
                                          printfn "pre-Cluster G = %f" preG
                                          newClusters
                                     else 
                                          x)
                        |> kmeanswapFn gainFn matrix depth
                        |> Seq.fold (fun (singles,best) i -> 
                                            let newNodes = 
                                                if singles=Map.empty then
                                                    i
                                                    |> Map.fold (fun state key nodes ->
                                                                        let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrix
                                                                        state 
                                                                        |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)
                                                else
                                                    i
                                                    |> Map.fold (fun state key nodes ->
                                                                        match (singles.TryFind key) with
                                                                        | None ->
                                                                            let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrix
                                                                            state 
                                                                            |> Map.add key (loop nodes (depth+1) dPredSum')
                                                                        | Some x ->
                                                                            state 
                                                                            |> Map.add key x                                                                   
                                                                    ) (Map.empty)
                                                                    
                                            let best' =
                                                if (confGainFn newNodes) > (confGainFn best) then  // compare configuration gains to get the best
                                                    newNodes  
                                                else 
                                                    best
                                            if (singles = Map.empty) then
                                                (newNodes, best')
                                            else
                                                (singles, best')
                                            
                                        ) (Map.empty,Map.empty) // here as state should be this singles (first) saved and optimal conf
                        |> snd

            
            {
            StepGain = stepGain; 
            ConfGain = (confGainFn children, children |> Map.toList |> List.map fst);
            GroupGain = max stepGain (confGainFn children);
            Member = nodeMembers;
            Children = 
                if  (mode=MM) || ((confGainFn children)>stepGain) then 
                    children;
                else 
                    Map.empty
            }
    
        loop rootGroup 0 0.   


    let matrixToPQ (sim: float [,]) =
        let n = sim.GetLength(1)
        Array.init (n) (fun id ->
            let p = MaxIndexPriorityQueue<float>(n) // n-i
            for j=0 to n-1 do 
                p.Insert j sim.[id,j]
            p
        )

    let createTreePQ gainFn kmeanswapFn clusterFn (weight: seq<float> option) (mode: Mode) (rootGroup: Types.Item array) = 

            let nRoot = rootGroup.Length

            let matrix = 
                rootGroup
                |> distMatrixWeightedOf distanceMatrixWeighted weight

            let pq = matrixToPQ matrix
                
            // calculation for one node    
            let rec loop (nodeMembers: Types.Item array) depth dPredSum =
        
                let idsParent = (groupIDFn nodeMembers)
                nodeMembers |> Array.iter (fun p -> pq.[p.ID].LeaveGroup idsParent)
                /// sum of max dist within a node 
                let dCurrSum = //dSumFn (groupIDFn nodeMembers) (groupIDFn nodeMembers) matrix
                    nodeMembers |> Array.sumBy (fun x -> pq.[x.ID].Top())

                /// to calc step from parent to current node
                let stepGain = (gainFn dCurrSum dPredSum nodeMembers.Length nRoot)

                let children = 
                    match mode with
                    |MM -> // MapMan with broken leaves
                        if nodeMembers.Length=1 then
                            Map.empty
                        else 
                            (breakGroup nodeMembers depth)
                            |> Map.map (fun key nodes -> 
                                            
                                            
                                            nodes |> Array.iter (fun p -> pq.[p.ID].LeaveGroup idsParent)
                                            let dPredSum' = nodes |> Array.sumBy (fun x -> pq.[x.ID].Top())
                                                        //dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrix
                                            //let pqArrayNew = pqArray |> Array.map (fun x -> x.DeepCopy())

                                            //nodes |> Array.iter (fun p -> pqArray.[p.ID].LeaveGroup (groupIDFn nodes)) 
                                            
                                            (loop nodes (depth+1) dPredSum'))

                    |SSN -> // with simplification overall
                        if      (nodeMembers.Length=1) 
                                || (nodeMembers.Length=0) 
                                || (String.contains "|" (String.Concat nodeMembers.[0].BinL))  then 
                            Map.empty
                        else 
                            (breakGroup nodeMembers depth) // get a map with original grouping as (subbin, Item)
                            |> (fun x -> 
                                        let singletons = x |> Map.filter (fun k pList -> pList.Length=1)
                                        if (singletons.Count>50) then
                                              let k = 50
                                              let list = singletons |> Map.toArray |> Array.map (snd) |> Array.concat |> List.ofArray
                                              let newClusters = clusterFn k weight list //mapOfSOM (snd (applySOM weight list 10 k)) list
                                              let oldClusters = x |> Map.filter (fun k pList -> pList.Length<>1)
                                              Map.fold (fun s k v -> Map.add k v s) newClusters oldClusters

                                         else 
                                              x)
                            |> kmeanswapFn gainFn matrix depth
                            |> Seq.fold (fun (singles,best) i -> 
                                                let newNodes = 
                                                    if singles=Map.empty then
                                                        i
                                                        |> Map.fold (fun state key nodes ->
                                                                            
                                                                nodes |> Array.iter (fun p -> pq.[p.ID].LeaveGroup idsParent)
                                                                let dPredSum' = nodes |> Array.sumBy (fun x -> pq.[x.ID].Top())
                                                                            //dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrix
                                                                
                                                                //let ids = (groupIDFn nodes)

                                                                //let pqArrayNew = pqArray |> Array.map (fun x -> x.DeepCopy())
                                                                //nodes |> Array.iter (fun p -> pqArray.[p.ID].LeaveGroup ids) 
                                                                            
                                                                state 
                                                                |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)
                                                    else
                                                        i
                                                        |> Map.fold (fun state key nodes ->
                                                            match (singles.TryFind key) with
                                                            | None ->

                                                                nodes |> Array.iter (fun p -> pq.[p.ID].LeaveGroup idsParent)
                                                                let dPredSum' = nodes |> Array.sumBy (fun x -> pq.[x.ID].Top())
                                                                                //dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrix
                                                                
                                                                //let ids = (groupIDFn nodes)

                                                                //let pqArrayNew = pqArray |> Array.map (fun x -> x.DeepCopy())
                                                                //nodes |> Array.iter (fun p -> pqArray.[p.ID].LeaveGroup ids) 

                                                                state 
                                                                |> Map.add key (loop nodes (depth+1) dPredSum')
                                                            | Some x ->
                                                                state 
                                                                |> Map.add key x                                                                   
                                                        ) (Map.empty)
                                                                    
                                                let best' =
                                                    if (confGainFn newNodes) > (confGainFn best) then  // compare configuration gains to get the best
                                                        newNodes  
                                                    else 
                                                        best
                                                if (singles = Map.empty) then
                                                    (newNodes, best')
                                                else
                                                    (singles, best')
                                            
                                            ) (Map.empty,Map.empty) // here as state should be this singles (first) saved and optimal conf
                            |> snd

            
                {
                StepGain = stepGain; 
                ConfGain = (confGainFn children, children |> Map.toList |> List.map fst);
                GroupGain = max stepGain (confGainFn children);
                Member = nodeMembers;
                Children = 
                    if  (mode=MM) || ((confGainFn children)>stepGain) then 
                        children;
                    else 
                        Map.empty
                }
    
            loop rootGroup 0 0.

    let createTreePQinKMS gainFn kmeanswapFn clusterFn (weight: seq<float> option) (mode: Mode) (rootGroup: Types.Item array) = 

            let nRoot = rootGroup.Length

            let matrix = 
                rootGroup
                |> distMatrixWeightedOf distanceMatrixWeighted weight

            let pq = matrixToPQ matrix
                
            // calculation for one node    
            let rec loop (nodeMembers: Types.Item array) depth dPredSum =
        
                let idsParent = (groupIDFn nodeMembers)
                nodeMembers |> Array.iter (fun p -> pq.[p.ID].LeaveGroup idsParent)
                /// sum of max dist within a node 
                let dCurrSum = //dSumFn (groupIDFn nodeMembers) (groupIDFn nodeMembers) matrix
                    nodeMembers |> Array.sumBy (fun x -> pq.[x.ID].Top())

                /// to calc step from parent to current node
                let stepGain = (gainFn dCurrSum dPredSum nodeMembers.Length nRoot)

                let children = 
                    match mode with
                    |MM -> // MapMan with broken leaves
                        if nodeMembers.Length=1 then
                            Map.empty
                        else 
                            (breakGroup nodeMembers depth)
                            |> Map.map (fun key nodes -> 
                                            nodes |> Array.iter (fun p -> pq.[p.ID].LeaveGroup idsParent)
                                            let dPredSum' = nodes |> Array.sumBy (fun x -> pq.[x.ID].Top())                                            
                                            (loop nodes (depth+1) dPredSum'))

                    |SSN -> // with simplification overall
                        if      (nodeMembers.Length=1) 
                                || (nodeMembers.Length=0) 
                                || (String.contains "|" (String.Concat nodeMembers.[0].BinL))  then 
                            Map.empty
                        else 
                            (breakGroup nodeMembers depth) // get a map with original grouping as (subbin, Item)
                            |> (fun x -> 
                                        let singletons = x |> Map.filter (fun k pList -> pList.Length=1)
                                        if (singletons.Count>50) then
                                              let k = 50
                                              let list = singletons |> Map.toArray |> Array.map (snd) |> Array.concat |> List.ofArray
                                              let newClusters = clusterFn k weight list
                                              let oldClusters = x |> Map.filter (fun k pList -> pList.Length<>1)
                                              Map.fold (fun s k v -> Map.add k v s) newClusters oldClusters

                                         else 
                                              x)
                            |> kmeanswapFn gainFn pq depth
                            |> Seq.fold (fun (singles,best) i -> 
                                                let newNodes = 
                                                    if singles=Map.empty then
                                                        i
                                                        |> Map.fold (fun state key nodes ->
                                                                            
                                                                nodes |> Array.iter (fun p -> pq.[p.ID].LeaveGroup idsParent)
                                                                let dPredSum' = nodes |> Array.sumBy (fun x -> pq.[x.ID].Top())
                                                                state 
                                                                |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)
                                                    else
                                                        i
                                                        |> Map.fold (fun state key nodes ->
                                                            match (singles.TryFind key) with
                                                            | None ->

                                                                nodes |> Array.iter (fun p -> pq.[p.ID].LeaveGroup idsParent)
                                                                let dPredSum' = nodes |> Array.sumBy (fun x -> pq.[x.ID].Top())
                                                                state 
                                                                |> Map.add key (loop nodes (depth+1) dPredSum')
                                                            | Some x ->
                                                                state 
                                                                |> Map.add key x                                                                   
                                                        ) (Map.empty)
                                                                    
                                                let best' =
                                                    if (confGainFn newNodes) > (confGainFn best) then  // compare configuration gains to get the best
                                                        newNodes  
                                                    else 
                                                        best
                                                if (singles = Map.empty) then
                                                    (newNodes, best')
                                                else
                                                    (singles, best')
                                            
                                            ) (Map.empty,Map.empty) // here as state should be this singles (first) saved and optimal conf
                            |> snd

            
                {
                StepGain = stepGain; 
                ConfGain = (confGainFn children, children |> Map.toList |> List.map fst);
                GroupGain = max stepGain (confGainFn children);
                Member = nodeMembers;
                Children = 
                    if  (mode=MM) || ((confGainFn children)>stepGain) then 
                        children;
                    else 
                        Map.empty
                }
    
            loop rootGroup 0 0.

    let getStepGainNodeSetnR setNR dCurrSum dPredSum numberCurr numberRoot =
        let nC = float numberCurr
        let nR = float setNR
        let deltaDist = dPredSum - dCurrSum
        let deltaSpec = -((nC/nR)*log2(nC/nR)+((nR-nC)/nR)*log2((nR-nC)/nR))
        if numberCurr=numberRoot then
            0.
        else
            deltaDist*deltaSpec

module PreClusterFunctions =
    
    
    /// call hierarchical clustering for singletons
    let clusterHier (k: int) weight (children: Types.Item list) =
        let clusters nClusters =
            children
            |> List.map (fun protein -> protein.dataL)
            |> ML.Unsupervised.HierarchicalClustering.generate (SSN.weightedEuclidean weight) (ML.Unsupervised.HierarchicalClustering.Linker.centroidLwLinker)
            |> ML.Unsupervised.HierarchicalClustering.cutHClust nClusters
            |> List.map (List.map (fun i -> children.[ML.Unsupervised.HierarchicalClustering.getClusterId i]))
        clusters k
        |> List.map (fun list ->
            let binName = list |> List.fold (fun acc x -> sprintf "%s|p%i" acc x.ID) "" 
            (binName,list |> Array.ofList))
        |> Map.ofList
    
    let rv2 (x: matrix) (y: matrix) =

        let xxt = x*x.Transpose
        let yyt = y*y.Transpose

        xxt |> Matrix.inplace_mapi (fun r c x -> if r = c then 0. else x)
        yyt |> Matrix.inplace_mapi (fun r c x -> if r = c then 0. else x)

        let num = (xxt*yyt).Diagonal |> Vector.sum
        let deno1 = xxt |> Matrix.map (fun x -> x**2.) |> Matrix.sum |> sqrt
        let deno2 = yyt |> Matrix.map (fun x -> x**2.) |> Matrix.sum |> sqrt

        num / (deno1 * deno2)

    let pairwiseCorrMax (x:matrix) (y:matrix) =
        let xN = x.Dimensions |> fst
        let yN = y.Dimensions |> fst
        let m = Array2D.create xN yN 0.
        for rowI in 0..(xN-1) do
            for colI in 0..(yN-1) do
                let tmp = SSN.weightedEuclidean None (x.Row rowI) (y.Row colI) 
                m.[rowI,colI] <- tmp
        m 
        |> Array2D.array2D_to_seq 
        |> Seq.max

    let pairwiseCorrAverage (x:matrix) (y:matrix) =
        let xN = x.Dimensions |> fst
        let yN = y.Dimensions |> fst
        let m = Array2D.create xN yN 0.
        for rowI in 0..(xN-1) do
            for colI in 0..(yN-1) do
                let tmp = SSN.weightedEuclidean None (x.Row rowI) (y.Row colI) 
                m.[rowI,colI] <- tmp
        m 
        |> Array2D.array2D_to_seq 
        |> Seq.average

    let clusterHierGroups (k: int) weight (children: Map<string,Types.Item []> ) =
        let clusters nClusters =
            let rvToMap (mapA: (string * Types.Item []))  (mapB: (string * Types.Item [])) = 
                let mA = mapA |> snd |> Array.map (fun protein -> protein.dataL) |> matrix
                let mB = mapB |> snd |> Array.map (fun protein -> protein.dataL) |> matrix
                pairwiseCorrAverage mA mB
            let children' = children |> Map.toArray
            children'
            |> ML.Unsupervised.HierarchicalClustering.generate rvToMap (ML.Unsupervised.HierarchicalClustering.Linker.completeLwLinker)
            |> ML.Unsupervised.HierarchicalClustering.cutHClust nClusters
            |> List.map (List.map (fun i -> children'.[ML.Unsupervised.HierarchicalClustering.getClusterId i]))
        clusters k
        |> List.map (fun list ->
            let binName = list |> List.fold (fun acc (bin,proteins) -> sprintf "%s|%s" acc bin) "" 
            let items = list |> List.map snd |> List.toArray |> Array.concat
            (binName,items))
        |> Map.ofList

    //// submodule about KMean for group clustering 

    let distMatrix (matrixA: string * matrix) (matrixB: string * matrix) =
        let mA = snd matrixA
        let mB = snd matrixB
        pairwiseCorrAverage mA mB

    let intiCgroups (input: (string*matrix) array) k : ((string*matrix) []) =

        let dmatrix = input |> Array.map (fun (bin,data) -> (data |> Matrix.enumerateColumnWise (fun x -> x |> Seq.median) )) |> MatrixTopLevelOperators.matrix
        
        //let cvmax = // find a feature with the biggest variance and return the (row Number, values of the feature), sorted by the values 
        //    dmatrix
        //    |> Matrix.Generic.enumerateColumnWise Seq.var
        //    |> Seq.zip (Matrix.Generic.enumerateColumnWise id dmatrix)
        //    |> Seq.maxBy snd
        //    |> fst
        //    |> Seq.mapi (fun rowI value -> (rowI,value)) 
        //    |> Seq.toArray 
        //    |> Array.sortBy snd

        let cvmax = // set the most important feature at 24H (with index 1) 
            dmatrix
            |> Matrix.Generic.enumerateColumnWise Seq.var
            |> Seq.zip (Matrix.Generic.enumerateColumnWise id dmatrix)
            |> Seq.item 1
            |> fst
            |> Seq.mapi (fun rowI value -> (rowI,value)) 
            |> Seq.toArray 
            |> Array.sortBy snd

        if cvmax.Length < k then failwithf "Number of data points must be at least %i" k        
        let chunkSize = cvmax.Length / k
        let midChunk  = chunkSize / 2
        [ for i=1 to k do
            let index = 
                match (chunkSize * i) with
                | x when x < cvmax.Length -> x - midChunk
                | x                       -> chunkSize * (i - 1) + ((cvmax.Length - chunkSize * (i - 1)) / 2)
            //printfn "Array.lenght = %i and index = %i" cvmax.Length (index-1)
            yield cvmax.[index-1] |> fst]
        |> Seq.map (fun rowI -> input.[rowI])
        |> Seq.toArray

    /// KKZ deterministic centroid initialization for kMean clustering
    let initCgroupsKKZ (data: (string*matrix) array) k =
        let centroid1 =
            data 
            |> Array.maxBy 
                (fun x -> 
                    x 
                    |> snd
                    |> Matrix.enumerateColumnWise (Seq.mean) 
                    |> fun xx ->  sqrt ( xx |> Seq.sumBy (fun i -> i*i))
                )
        let LeaveData d c =
            d |> Array.removeIndex (Array.FindIndex<string*matrix>(d, fun x -> x=c))

        let rec loop dataRest kRest centroids =
            if kRest=1 then   
                centroids
            else    
                let newC = 
                    dataRest 
                    |> Array.map (fun (s,p) -> 
                        (s,p), centroids |> List.map (fun (sc,c) -> pairwiseCorrAverage (p)  (c)) |> List.min )
                    |> Array.maxBy snd 
                    |> fst
                loop (LeaveData dataRest newC) (kRest-1) (newC::centroids)

        loop (LeaveData data centroid1) k [centroid1]
        |> List.toArray

    // Recompute Centroid as average of given sample (for kmeans)
    let updateCentroid (current: string * matrix) (sample: (string * matrix) []) = // rewrite it in matrix!
        let size = sample.Length
        match size with
        | 0 -> current
        | _ ->
            ("", 
                sample
                |> Array.map (fun (_,x) -> x.ToArray2D() |> Array2D.toJaggedArray)
                |> Array.concat
                |> fun p -> MatrixTopLevelOperators.matrix p)

    //    /// Convert centroids into an initial scheme for K-Mean-Swap
    let centroidsToScheme (input: Item [] []) (centroid: matrix []) (scheme: int []) : ((string*(Item [])) []) =
        let sorted = input |> Array.indexed |>  Array.sortByDescending ( fun (i,x) -> Array.length x)
        let rec loop itemsRest groupID =
            [|if groupID=scheme.Length then   
                let binName = itemsRest |> Array.concat |> Array.map (fun p -> p.ID) |> Array.sort |> fun i -> String.Join("|",i)
                yield (binName,itemsRest |> Array.concat)
            else    
                let (itemsCluster,itemsNew) = 
                    itemsRest 
                    |> Array.sortByDescending 
                        (fun i -> pairwiseCorrAverage centroid.[groupID] (i |> Array.map (fun x -> x.dataL) |> MatrixTopLevelOperators.matrix))
                    |> Array.splitAt scheme.[groupID]
                let binName = itemsCluster |> Array.concat |> Array.map (fun p -> p.ID) |> Array.sort |> fun i -> String.Join("|",i)
                yield (binName,itemsCluster |> Array.concat)
                yield! loop itemsNew (groupID+1)
            |]
        loop (input) 0

    let clustersToCentroidMatrix (clusters: Item [] [] []) : matrix [] =
        clusters
        |> Array.map (fun cluster -> 
            cluster
            |> Array.map (fun group ->
                group |> Array.map (fun x -> x.dataL)
                )
            |> Array.concat
            |> MatrixTopLevelOperators.matrix
            )

    let centroidsToMatrix (centroids: Item [] list) : matrix [] =
        centroids
        |> List.toArray
        |> Array.map (fun centroid ->
            centroid |> Array.map (fun x -> x.dataL) |> MatrixTopLevelOperators.matrix
            )

    let kkzSingle (data: Item []) k =
        let centroid1 =
            data |> Array.maxBy (fun x -> sqrt ( x.dataL |> Array.sumBy (fun i -> i*i)))
        let LeaveData d c =
            d |> Array.removeIndex (Array.FindIndex<Item>(d, fun x -> x=c))
        let rec loop dataRest kRest centroids =
            if kRest=1 then   
                centroids
            else    
                let newC = 
                    dataRest 
                    |> Array.map  (fun p -> p, centroids |> List.map (fun c -> SSN.weightedEuclidean None p.dataL c.dataL) |> List.min )
                    |> Array.maxBy snd 
                    |> fst
                loop (LeaveData dataRest newC) (kRest-1) (newC::centroids)
        loop (LeaveData data centroid1) k [centroid1]

    /// give a list of centroids as k the most distant elements of the dataset     
    let kkz (data: Item [] []) k =
        let centroid1 =
            data 
            |> Array.maxBy 
                (fun x -> 
                    x 
                    |> Array.map (fun x -> x.dataL) 
                    |> JaggedArray.transpose 
                    |> Array.map (Array.average) 
                    |> fun xx ->  sqrt ( xx |> Array.sumBy (fun i -> i*i))
                )
        let LeaveData d c =
            d |> Array.removeIndex (Array.FindIndex<Item []>(d, fun x -> x=c))
        let toMatrix =
            Array.map (fun i -> i.dataL)
            >> MatrixTopLevelOperators.matrix

        let rec loop dataRest kRest centroids =
            if kRest=1 then   
                centroids
            else    
                let newC = 
                    dataRest 
                    |> Array.map (fun p -> 
                        p, centroids |> List.map (fun c -> pairwiseCorrAverage (toMatrix p)  (toMatrix c)) |> List.min )
                    |> Array.maxBy snd 
                    |> fst
                loop (LeaveData dataRest newC) (kRest-1) (newC::centroids)
        loop (LeaveData data centroid1) k [centroid1]

    let varPart_Single (data: Item []) k =
        let sse (cluster: Item []) =
            let centroid = cluster |> Array.map (fun x -> x.dataL) |> MatrixTopLevelOperators.matrix
            let dist (a: Item) (b: matrix) =
                pairwiseCorrAverage ([a.dataL] |> MatrixTopLevelOperators.matrix) b
            cluster
            |> Array.map 
                (fun i -> (dist i centroid)*(dist i centroid))
            |> Array.average
        let split (cluster: Item []) =
            let featureN = 
                cluster 
                |> Array.map (fun x -> x.dataL) 
                |> MatrixTopLevelOperators.matrix 
                |> Matrix.Generic.enumerateColumnWise Seq.var
                |> Seq.mapi (fun id x -> (id,x))
                |> Seq.maxBy snd
                |> fst
            let featureMean =
                cluster
                |> Array.map (fun x -> x.dataL.[featureN])
                |> Array.average
            cluster
            |> Array.partition (fun x -> x.dataL.[featureN]>featureMean)
        let pq = MaxIndexPriorityQueue<float>(k*2-1) // put a cluster there or SSE of a cluster or SSE*cluster????
        let clusters: Item [] [] = Array.create (k*2-1) [||]
        pq.Insert 0 (sse data)
        clusters.[0] <- data
        [|1 .. (k-1)|] 
        |> Array.iter (fun ik -> 
            let loosest = clusters.[pq.HeapItemIndex 1]
            let newCl = split loosest
            pq.Pop() |> ignore
            pq.Insert (2*ik) (sse (fst newCl))
            pq.Insert (2*ik-1) (sse (snd newCl))
            clusters.[2*ik] <- (fst newCl)
            clusters.[2*ik-1] <- (snd newCl)
            )
        [|1 .. k|]
        |> Array.map (fun x -> clusters.[pq.HeapItemIndex x])


    let varPart (data: Item [] []) k =
        let sse (cluster: Item [] []) =
            let centroid = cluster |> Array.map (fun x -> x |> Array.map (fun i -> i.dataL)) |> Array.concat |> MatrixTopLevelOperators.matrix
            let dist (a: Item []) (b: matrix) =
                pairwiseCorrAverage (a |> Array.map (fun i -> i.dataL) |> MatrixTopLevelOperators.matrix) b
            cluster
            |> Array.map 
                (fun i -> (dist i centroid)*(dist i centroid))
            |> Array.average
        let split (cluster: Item [] []) =
            let featureN = 
                cluster 
                |> Array.map (fun x -> x |> Array.map (fun i -> i.dataL) |> JaggedArray.transpose |> Array.map Array.average )
                |> MatrixTopLevelOperators.matrix 
                |> Matrix.Generic.enumerateColumnWise Seq.var
                |> Seq.mapi (fun id x -> (id,x))
                |> Seq.maxBy snd
                |> fst
            let featureMean =
                cluster
                |> Array.map (fun x -> x |> Array.map (fun i -> i.dataL.[featureN]) |> Array.average)
                |> Array.average
            cluster
            |> Array.partition (fun x -> (x |> Array.map (fun i -> i.dataL.[featureN]) |> Array.average)>featureMean)
        let pq = MaxIndexPriorityQueue<float>(k*2-1) // put a cluster there or SSE of a cluster or SSE*cluster????
        let clusters: Item [] [] [] = Array.create (k*2-1) [||]
        pq.Insert 0 (sse data)
        clusters.[0] <- data
        [|1 .. (k-1)|] 
        |> Array.iter (fun ik -> 
            let loosest = clusters.[pq.HeapItemIndex 1]
            let newCl = split loosest
            pq.Pop() |> ignore
            pq.Insert (2*ik) (sse (fst newCl))
            pq.Insert (2*ik-1) (sse (snd newCl))
            clusters.[2*ik] <- (fst newCl)
            clusters.[2*ik-1] <- (snd newCl)
            )
        [|1 .. k|]
        |> Array.map (fun x -> clusters.[pq.HeapItemIndex x])

    let kmeanGroups (k: int) weight (children: Map<string,Types.Item []> ) : Map<string,Types.Item []> =

        let data = children |> Map.toArray |> Array.map (fun (s,ar) -> (s,ar |> Array.map (fun p -> p.dataL) |> MatrixTopLevelOperators.matrix))

        let clusters = 
            let c1 = ML.Unsupervised.IterativeClustering.compute distMatrix (intiCgroups) updateCentroid data k
            let x1 = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c1
            [|1 .. 20|]
            |> Array.fold (fun (disp,best) x -> 
                let c = ML.Unsupervised.IterativeClustering.compute distMatrix (intiCgroups) updateCentroid data k
                let x = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c
                if x<disp then
                    (x,c)
                else
                    (disp,best) ) (x1,c1)
            |> snd

        data
        |> Array.map (fun list -> (clusters.Classifier list |> fst),list )
        |> Array.groupBy (fst)
        |> Array.map (fun (cID,list) ->
            let binName = list |> Array.map (fun (cID,(bin,p)) -> bin) |> Array.sort |> fun i -> String.Join("|",i) 
            let items = list |> Array.map (fun (cID,(bin,p)) -> children |> Map.find bin) |> Array.concat
            (binName,items))
        |> Map.ofArray
    
    let centroidFactory (input: (string*matrix) array) (k: int) : ((string*matrix) []) =
        let r = new System.Random() 
        ML.Unsupervised.IterativeClustering.randomCentroids r input k

    let kmeanGroupsRandom (k: int) weight (children: Map<string,Types.Item []> ) : Map<string,Types.Item []> =

        let data = children |> Map.toArray |> Array.map (fun (s,ar) -> (s,ar |> Array.map (fun p -> p.dataL) |> MatrixTopLevelOperators.matrix))

        let clusters = 
            let c1 = ML.Unsupervised.IterativeClustering.compute distMatrix (centroidFactory) updateCentroid data k
            let x1 = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c1
            [|1 .. 20|]
            |> Array.fold (fun (disp,best) x -> 
                let c = ML.Unsupervised.IterativeClustering.compute distMatrix (centroidFactory) updateCentroid data k
                let x = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c
                if x<disp then
                    (x,c)
                else
                    (disp,best) ) (x1,c1)
            |> snd

        data
        |> Array.map (fun list -> (clusters.Classifier list |> fst),list )
        |> Array.groupBy (fst)
        |> Array.map (fun (cID,list) ->
            let binName = list |> Array.map (fun (cID,(bin,p)) -> bin) |> Array.sort |> fun i -> String.Join("|",i) 
            let items = list |> Array.map (fun (cID,(bin,p)) -> children |> Map.find bin) |> Array.concat
            (binName,items))
        |> Map.ofArray



module KMeanSwapFunctions =
   
    let gainstepcheck f (matrix: float [,]) (conf: int array array)  =
        let parentGroup = conf |> Array.concat
        let gFn (current: int array) = SSN.getStepGainFn f current parentGroup parentGroup.Length matrix
        conf 
        |> Array.fold (fun acc x -> acc + (gFn x)) 0.

    let gainstepcheckPQ f (pq) (conf: int array array)  =
        let parentGroup = conf |> Array.concat
        let gFn (current: int array) = SSN.getStepGainFnPQ f current parentGroup parentGroup.Length pq
        conf 
        |> Array.fold (fun acc x -> acc + (gFn x)) 0.

///////////////// old function: arrays, fixed initial values, no saving the best configuration

    let kMeanSwapOld f matrixItems depth (initialGrouping: (string*(Types.Item [])) []) : (Map<string,Types.Item []>) seq = 
    
        let parentGroup = initialGrouping |> Array.map (fun (k,p) -> p) |> Array.concat

        //let mutable iterN : int = 0

        let schemes =  seq [2..initialGrouping.Length] |> Seq.map (SSN.schemeGenerator initialGrouping.Length) |> Seq.concat

        let gainDiff listA listB : float =
            let gFn (current: Types.Item []) = SSN.getStepGainFn f (SSN.groupIDFn current) (SSN.groupIDFn parentGroup) parentGroup.Length matrixItems
            (gFn (Array.append listA listB))

        let matrix =     
                let data =
                    initialGrouping 
                    |> Array.map snd
                let m = Array2D.zeroCreate (data.Length) (data.Length)
                for rowI in 0..data.Length-1 do
                    for colI in 0..rowI do
                        let tmp = if rowI=colI then nan else gainDiff data.[rowI] data.[colI]
                        m.[colI,rowI] <- tmp
                        m.[rowI,colI] <- tmp
                m

        let processInt (intListList: int [] []) =
            intListList
            |> Array.map (fun i ->
                                    if i.Length=1 then
                                        initialGrouping.[i.[0]]
                                    else
                                        Array.fold (fun (key,gr) ii ->
                                        let (binKey, value) = initialGrouping.[ii]
                                        ((sprintf "%s|%s" key binKey), Array.append value gr )
                                        ) ("",[||]) i)
            |> Array.map (fun (newBin,protA) ->
                (newBin, protA |> Array.map (fun prot ->
                        {prot with BinL = Array.append prot.BinL.[0..(depth)] [|newBin|]})))  
            |> Map.ofArray
    
        let rec swap (cluster: int []) (currentConf: int [] []) n = /// K-Mean-Swap for matrices
            
                let outside = 
                    currentConf
                    |> Array.concat
                    |> Array.except cluster
                let (idOut, idIn, minDist) =
                    outside
                    |> Array.map (fun i -> Array.map (fun j -> (i, j, matrix.[i,j])) cluster)
                    |> Array.concat
                    |> Array.maxBy (fun (outsideItem, insideItem, distance) -> distance) // find the closest item from outside, 
                let (farItem, closeItem, maxDist) =
                    cluster
                    |> Array.filter (fun iIn -> iIn<>idIn)
                    |> Array.map (fun i -> (i, idIn, matrix.[i,idIn])) 
                    |> Array.minBy (fun (farItem, closeItem, distance) -> distance) // find the farthest to the inside item within a cluster

                if (n >= 0)  then 
                    if (maxDist) < (minDist) then
                        let idInsideItem = cluster |> Array.findIndex (fun id -> id = farItem) 
                        let idOutsideCluster = currentConf |>  Array.findIndex (fun i -> i |> Array.exists (fun (ii) -> ii = idOut))
                        let idOutsideItem = Array.findIndex (fun id -> id = idOut) currentConf.[idOutsideCluster]
                        let temp = cluster.[idInsideItem]
                        cluster.[idInsideItem] <- currentConf.[idOutsideCluster].[idOutsideItem]
                        currentConf.[idOutsideCluster].[idOutsideItem] <- temp
                        //iterN <- iterN + 1
                        swap cluster currentConf (n-1) 
                //printfn "quit swap loop"

        let fillSchemes (scheme: int []) : Map<string, Types.Item []> =
    
            let mutable configurationIDs = 
                scheme
                |> Array.mapFold (fun before n -> ([|before .. (before+n-1)|],(before+n))) 0
                |> fst

            configurationIDs
            |> Array.iter (fun cluster -> if (cluster.Length>1) then swap cluster configurationIDs cluster.Length ) // not saving the best
        
            let res =
                configurationIDs
                |> processInt
            res

        //let r =                                   if want to see iteration count
        //    schemes
        //    |> Seq.rev
        //    |> List.ofSeq
        //    |> List.map fillSchemes
        //    |> Seq.ofList
        //printfn "iterations done: %i" iterN
        //r

        schemes
        |> Seq.rev
        |> Seq.map fillSchemes
    
    let kmeanSwapShuffleOld power f matrixItems depth (map : Map<string,Types.Item []>) =
        let items = map |> Map.toArray
        //kMeanSwapOld f matrixItems depth items
        //|> Seq.toArray
        [|for a in [1 .. power] do yield
                                        items
                                        |> Array.shuffleFisherYates // fn how to initialize first groups
                                        |> (fun x -> 
                                                printfn "randomWalk %i" a
                                                kMeanSwapOld f matrixItems depth (x))
                                        |> Seq.toArray|]
        |> Array.stackVertical
        |> Array2D.toJaggedArray
        |> Array.map (fun a -> a |> Array.maxBy (fun m -> m |> Map.toArray |> Array.map (snd >> SSN.groupIDFn) |> gainstepcheck f matrixItems))

///////////////////////// Best Random Walk - arrays, random initial values, save best from configuration in a random walk

    let kMeanSwapNew f matrixItems depth (initialGrouping: (string*(Types.Item array)) []) : (Map<string,Types.Item array>) seq =
    
        let parentGroup = initialGrouping |> Array.map (fun (k,p) -> p) |> Array.concat

        let mutable iterN : int = 0

        let schemes =  seq [2..initialGrouping.Length] |> Seq.map (SSN.schemeGenerator initialGrouping.Length) |> Seq.concat

        let gainDiff listA listB : float =
            let gFn (current: Types.Item array) = SSN.getStepGainFn f (SSN.groupIDFn current) (SSN.groupIDFn parentGroup) parentGroup.Length matrixItems
            (gFn (Array.append listA listB))

        let matrix =     
            let data =
                initialGrouping
                |> Array.map snd
            let m = Array2D.zeroCreate (data.Length) (data.Length)
            for rowI in 0..data.Length-1 do
                for colI in 0..rowI do
                    let tmp = if rowI=colI then nan else gainDiff data.[rowI] data.[colI]
                    m.[colI,rowI] <- tmp
                    m.[rowI,colI] <- tmp
            m

        let processInt (intListList: int [] []) =
            intListList
            |> Array.map (fun i ->
                                    if i.Length=1 then
                                        initialGrouping.[i.[0]]
                                    else
                                        Array.fold (fun (key,gr) ii ->
                                        let (binKey, value) = initialGrouping.[ii]
                                        ((sprintf "%s|%s" key binKey), Array.append value gr )
                                        ) ("",[||]) i)
            |> Array.map (fun (newBin,protA) ->
                (newBin, protA |> Array.map (fun prot ->
                        {prot with BinL = Array.append prot.BinL.[0..(depth)] [|newBin|]}))) // updating the bin labels            
            |> Map.ofArray

        let rec swap (matrix: float [,]) (clusterN: int) (currentConf: int [] []) n downIterLast bestSoFar = /// K-Mean-Swap for matrices
            
            let config = currentConf |> Array.map (Array.map (fun i -> initialGrouping.[i] |> snd |> SSN.groupIDFn) >> Array.concat)
            let evalGain = gainstepcheck f matrixItems config //////// optimize the recalculation!

            let (downIter, best) = 
                if (evalGain < (fst bestSoFar)) then
                    (downIterLast+1, bestSoFar)
                else
                    (downIterLast, (evalGain,currentConf))
        
            //accResults <- (evalGain,temp)::accResults

            let cluster = currentConf.[clusterN]

            let outside =
                currentConf
                |> Array.concat
                |> Array.except cluster
            let (idOut, idIn, minDist) =
                outside
                |> Array.map (fun iOut -> Array.map (fun iIn -> (iOut, iIn, matrix.[iOut,iIn])) cluster)
                |> Array.concat
                |> Array.maxBy (fun (_, _, distance) -> distance) // find the closest item from outside,
            let (farItem, _, maxDist) =
                cluster
                |> Array.filter (fun iIn -> iIn<>idIn)
                |> Array.map (fun iIn -> (iIn, idIn, matrix.[iIn,idIn]))
                |> Array.minBy (fun (_, _, distance) -> distance) // find the farthest to the inside item within a cluster
            
            if (n >= 0) // size of the clsuter as a limit for amount of iterations
                && (downIter<2) // no more than 2 down-iterations pro cluster
                && (maxDist) < (minDist) // there is a gain in swapping
                then
                    let idInsideItem = cluster |> Array.findIndex (fun id -> id = farItem)
                    let idOutsideCluster = currentConf |>  Array.findIndex (fun i -> i |> Array.exists (fun (ii) -> ii = idOut))
                    let idOutsideItem = Array.findIndex (fun id -> id = idOut) currentConf.[idOutsideCluster]
                    let temp = cluster.[idInsideItem]
                    currentConf.[clusterN].[idInsideItem] <- currentConf.[idOutsideCluster].[idOutsideItem]
                    currentConf.[idOutsideCluster].[idOutsideItem] <- temp
                    //printfn "config: %A" currentConf
                    iterN <- iterN + 1
                    swap matrix clusterN currentConf (n-1) downIter best
                else best

        let fillSchemes (scheme: int []) : Map<string, Types.Item array> =
        
            let configurationIDs =
                scheme
                |> Array.mapFold (fun before n -> ([|before .. (before+n-1)|],(before+n))) 0
                |> fst
        
            // accResults <- []

            let result =
                [|0 .. (scheme.Length-1)|]   
                |> Array.fold (fun (gainBest, bestSoFar) clusterN -> 
                    if (scheme.[clusterN]>1) then 
                        swap matrix clusterN bestSoFar scheme.[clusterN] 0 (0.,[||]) 
                    else 
                        (0.,bestSoFar)) (0., configurationIDs)
                |> snd
            
            result |> processInt
    
        let r = 
            schemes
            |> Seq.rev
            |> List.ofSeq
            |> List.map fillSchemes
            |> Seq.ofList

        printfn "iterations done: %i" iterN

        r

    let kmeanSwapShuffleNew power f matrixItems depth (map : Map<string,Types.Item array>) =
        let items = map |> Map.toArray
        [|for a in [1 .. power] do yield
                                        items
                                        //|> Array.shuffleFisherYates 
                                        |> (fun x -> 
                                                printfn "randomWalk %i" a
                                                kMeanSwapNew f matrixItems depth x) //|> Map.ofArray))
                                        |> Seq.toArray|]
        |> Array.stackVertical
        |> Array2D.toJaggedArray
        |> Array.map (fun a -> a |> Array.maxBy (fun m -> m |> Map.toArray |> Array.map (snd >> SSN.groupIDFn) |> gainstepcheck f matrixItems))
    
////////////////////////    Gain Calculation by components - arrays, best from several random walk, gain calculations carried during kMean swapping

    type GainComponents = {
        ParentDiss : float
        CurrentDiss : float
        }

    let gainValueFn nSet (gainA: GainComponents [] []) =
        gainA
        |> Array.fold (fun acc group ->
            let nC = float group.Length
            let nR = float nSet
            let deltaSpec = -((nC/nR)*log2(nC/nR)+((nR-nC)/nR)*log2((nR-nC)/nR))
            acc + (group |> Array.sumBy (fun x -> x.ParentDiss-x.CurrentDiss))*deltaSpec
            ) 0.

    /// save gain calculation by saving some gain calculation components intact
    let kMeanSwapSaveGain nSet f matrixSingletons depth (initialGrouping: (string*(Types.Item array)) []) : (Map<string,Types.Item array>) seq =
    
        let parentGroup = initialGrouping |> Array.map (fun (k,p) -> p) |> Array.concat

        let mutable iterN : int = 0

        let schemes =  seq [2..initialGrouping.Length] |> Seq.map (SSN.schemeGenerator initialGrouping.Length) |> Seq.concat

        let gainDiff listA listB : float =
            let gFn (current: Types.Item array) = 
                SSN.getStepGainFn f (SSN.groupIDFn current) (SSN.groupIDFn parentGroup) parentGroup.Length matrixSingletons
            (gFn (Array.append listA listB))

        let matrix =     
            let data =
                initialGrouping
                |> Array.map snd
            let m = Array2D.zeroCreate (data.Length) (data.Length)
            for rowI in 0..data.Length-1 do
                for colI in 0..rowI do
                    let tmp = if rowI=colI then nan else gainDiff data.[rowI] data.[colI]
                    m.[colI,rowI] <- tmp
                    m.[rowI,colI] <- tmp
            m

        let processInt (intListList: int [] []) =
            intListList
            |> Array.map (fun i ->
                                    if i.Length=1 then
                                        initialGrouping.[i.[0]]
                                    else
                                        Array.fold (fun (key,gr) ii ->
                                        let (binKey, value) = initialGrouping.[ii]
                                        ((sprintf "%s|%s" key binKey), Array.append value gr )
                                        ) ("",[||]) i)
            |> Array.map (fun (newBin,protA) ->
                (newBin, protA |> Array.map (fun prot ->
                        {prot with BinL = Array.append prot.BinL.[0..(depth)] [|newBin|]}))) // updating the bin labels            
            |> Map.ofArray

        let rec swap (matrix: float [,]) (clusterN: int) (currentConf: int [] []) n downIterLast bestSoFar (lastGain: GainComponents [] []) (swappedGroups: int []) (swappedElements: int []) = /// K-Mean-Swap for matrices
            
            let newGain =
                if swappedGroups=[||] then
                    lastGain
                else
                    lastGain
                    |> Array.mapi (fun idGroup group ->
                        if (idGroup = swappedGroups.[0]) then
                            group
                            |> Array.mapi (fun idElement element ->

                                let currentGroupIDs = currentConf.[idGroup] |> Array.map (fun i -> initialGrouping.[i] |> snd) |> Array.concat  |> SSN.groupIDFn
                                let currentElement = initialGrouping.[currentConf.[idGroup].[idElement]] |> snd |> SSN.groupIDFn

                                if idElement = swappedElements.[0] then
                                    {ParentDiss = SSN.dSumFn currentElement (parentGroup |> SSN.groupIDFn) matrixSingletons;
                                    CurrentDiss = SSN.dSumFn currentElement currentGroupIDs matrixSingletons}
                                else
                                    {element with CurrentDiss = SSN.dSumFn currentElement currentGroupIDs matrixSingletons})
                        elif (idGroup = swappedGroups.[1]) then
                            group
                            |> Array.mapi (fun idElement element ->
                        
                                let currentGroupIDs = currentConf.[idGroup] |> Array.map (fun i -> initialGrouping.[i] |> snd) |> Array.concat  |> SSN.groupIDFn
                                let currentElement = initialGrouping.[currentConf.[idGroup].[idElement]] |> snd |> SSN.groupIDFn

                                if idElement = swappedElements.[1] then
                                    {ParentDiss = SSN.dSumFn currentElement (parentGroup |> SSN.groupIDFn) matrixSingletons;
                                     CurrentDiss = SSN.dSumFn currentElement currentGroupIDs matrixSingletons}
                                else
                                    {element with CurrentDiss = SSN.dSumFn currentElement currentGroupIDs matrixSingletons})
                        else
                            group
                        )

            let newGainValue =
                gainValueFn nSet newGain

            let (downIter, best) =
                if (newGainValue < (fst bestSoFar)) then
                    (downIterLast+1, bestSoFar)
                else
                    (downIterLast, (newGainValue,currentConf))
        
            let cluster = currentConf.[clusterN]

            let outside =
                currentConf
                |> Array.concat
                |> Array.except cluster
            let (idOut, idIn, minDist) =
                outside
                |> Array.map (fun iOut -> Array.map (fun iIn -> (iOut, iIn, matrix.[iOut,iIn])) cluster)
                |> Array.concat
                |> Array.maxBy (fun (_, _, distance) -> distance) // find the closest item from outside,
            let (farItem, _, maxDist) =
                cluster
                |> Array.filter (fun iIn -> iIn<>idIn)
                |> Array.map (fun iIn -> (iIn, idIn, matrix.[iIn,idIn]))
                |> Array.minBy (fun (_, _, distance) -> distance) // find the farthest to the inside item within a cluster
            
            if (n >= 0) // size of the cluster as a limit for amount of iterations
                && (downIter<2) // no more than 2 down-iterations pro cluster
                && (maxDist) < (minDist) // there is a gain in swapping
                then
                    let idInsideItem = cluster |> Array.findIndex (fun id -> id = farItem)
                    let idOutsideCluster = currentConf |>  Array.findIndex (fun i -> i |> Array.exists (fun (ii) -> ii = idOut))
                    let idOutsideItem = Array.findIndex (fun id -> id = idOut) currentConf.[idOutsideCluster]
                    let temp = cluster.[idInsideItem]
                    currentConf.[clusterN].[idInsideItem] <- currentConf.[idOutsideCluster].[idOutsideItem]
                    currentConf.[idOutsideCluster].[idOutsideItem] <- temp
                    iterN <- iterN + 1
                    swap matrix clusterN currentConf (n-1) downIter best newGain [|clusterN;idOutsideCluster|] [|idInsideItem;idOutsideItem|]
                else best

        let fillSchemes (scheme: int []) : Map<string, Types.Item array> =
        
            let configurationIDs =
                scheme
                |> Array.mapFold (fun before n -> ([|before .. (before+n-1)|],(before+n))) 0
                |> fst
        
            let result =
                [|0 .. (scheme.Length-1)|]
                |> Array.fold (fun (gainBest, bestSoFar) clusterN ->
                    if (scheme.[clusterN]>1) then
                    
                        let gainComp =
                            configurationIDs
                            |> Array.map (fun idGroup ->
                                Array.map (fun idElement ->
                                    let currentGroupIDs = idGroup |> Array.map (fun i -> initialGrouping.[i] |> snd) |> Array.concat  |> SSN.groupIDFn
                                    let currentElement = initialGrouping.[idElement] |> snd |> SSN.groupIDFn
                                    {ParentDiss = SSN.dSumFn currentElement (parentGroup |> SSN.groupIDFn) matrixSingletons;
                                    CurrentDiss = SSN.dSumFn currentElement currentGroupIDs matrixSingletons}) idGroup )

                        swap matrix clusterN bestSoFar scheme.[clusterN] 0 (0.,[||]) gainComp [||] [||]
                    else
                        (gainBest,bestSoFar)) (0., configurationIDs)
                |> snd
            
            result |> processInt
    
        let r =
            schemes
            |> Seq.rev
            |> List.ofSeq
            |> List.map fillSchemes
            |> Seq.ofList

        //printfn "iterations done: %i" iterN

        r

    let kmeanSwapShuffle setN power f matrixItems depth (map : Map<string,Types.Item array>) =
        let items = map |> Map.toArray
        [|for a in [1 .. power] do yield
                                        items
                                        |> Array.shuffleFisherYates
                                        |> (fun x ->
                                                //printfn "randomWalk %i" a
                                                kMeanSwapSaveGain setN f matrixItems depth x )
                                        |> Seq.toArray|]
        |> Array.stackVertical
        |> Array2D.toJaggedArray
        |> Array.map (fun a -> a |> Array.maxBy (fun m -> m |> Map.toArray |> Array.map (snd >> SSN.groupIDFn) |> gainstepcheck f matrixItems))


    /// Save gain calculation and use Priority Queue instead of singletonMatrix
    let kMeanSwapSaveGainPQ nSet f pq depth (initialGrouping: (string*(Types.Item array)) []) : (Map<string,Types.Item array>) seq =
    
        let getStepGainFnPQ fn itemsChild itemsParent (nRoot: int) (pq: MaxIndexPriorityQueue<float> []) : float =
        
            let dPredSum = itemsChild |> Array.sumBy (fun x -> 
                                                            pq.[x].LeaveGroup itemsParent
                                                            pq.[x].Top())
            let dCurrSum = itemsChild |> Array.sumBy (fun x -> 
                                                            pq.[x].LeaveGroup itemsChild
                                                            pq.[x].Top())
            fn dCurrSum dPredSum itemsChild.Length nRoot 

        let parentGroup = initialGrouping |> Array.map (fun (k,p) -> p) |> Array.concat
        let idsParent = SSN.groupIDFn parentGroup

        let mutable iterN : int = 0

        let schemes =  seq [2..initialGrouping.Length] |> Seq.map (SSN.schemeGenerator initialGrouping.Length) |> Seq.concat

        let gainDiff listA listB : float =
            let gFn (current: Types.Item array) = getStepGainFnPQ f (SSN.groupIDFn current) idsParent parentGroup.Length pq
            (gFn (Array.append listA listB))

        let matrixGainsForGroups =     
            let data =
                initialGrouping
                |> Array.map snd
            let m = Array2D.zeroCreate (data.Length) (data.Length)
            for rowI in 0..data.Length-1 do
                for colI in 0..rowI do
                    let tmp = if rowI=colI then nan else gainDiff data.[rowI] data.[colI]
                    m.[colI,rowI] <- tmp
                    m.[rowI,colI] <- tmp
            m

        let processInt (intListList: int [] []) =
            intListList
            |> Array.map (fun i ->
                                    if i.Length=1 then
                                        initialGrouping.[i.[0]]
                                    else
                                        Array.fold (fun (key,gr) ii ->
                                        let (binKey, value) = initialGrouping.[ii]
                                        ((sprintf "%s|%s" key binKey), Array.append value gr )
                                        ) ("",[||]) i)
            |> Array.map (fun (newBin,protA) ->
                (newBin, protA |> Array.map (fun prot ->
                        {prot with BinL = Array.append prot.BinL.[0..(depth)] [|newBin|]}))) // updating the bin labels            
            |> Map.ofArray

        let rec swap (matrix: float [,]) (clusterN: int) (currentConf: int [] []) n downIterLast bestSoFar (lastGain: GainComponents [] []) (swappedGroups: int []) (swappedElements: int []) = /// K-Mean-Swap for matrices

            let newGain =
                if swappedGroups=[||] then
                    lastGain
                else
                    lastGain
                    |> Array.mapi (fun idGroup group ->
                        if (idGroup = swappedGroups.[0]) then
                            group
                            |> Array.mapi (fun idElement element ->

                                let currentGroupIDs = currentConf.[idGroup] |> Array.map (fun i -> initialGrouping.[i] |> snd) |> Array.concat  |> SSN.groupIDFn
                                let currentElement = initialGrouping.[currentConf.[idGroup].[idElement]] |> snd |> SSN.groupIDFn

                                if idElement = swappedElements.[0] then
                                    let dPredSum = currentElement |> Array.sumBy (fun x -> 
                                                                                    pq.[x].LeaveGroup idsParent
                                                                                    pq.[x].Top())
                                    let dCurrSum = currentElement |> Array.sumBy (fun x -> 
                                                                                    pq.[x].LeaveGroup currentGroupIDs
                                                                                    pq.[x].Top())
                                    {ParentDiss = dPredSum;
                                    CurrentDiss = dCurrSum}
                                else
                                    let dCurrSum = currentElement |> Array.sumBy (fun x -> 
                                                                                    pq.[x].LeaveGroup currentGroupIDs
                                                                                    pq.[x].Top())
                                    {element with CurrentDiss = dCurrSum})
                        elif (idGroup = swappedGroups.[1]) then
                            group
                            |> Array.mapi (fun idElement element ->
                        
                                let currentGroupIDs = currentConf.[idGroup] |> Array.map (fun i -> initialGrouping.[i] |> snd) |> Array.concat  |> SSN.groupIDFn
                                let currentElement = initialGrouping.[currentConf.[idGroup].[idElement]] |> snd |> SSN.groupIDFn

                                if idElement = swappedElements.[1] then
                                    let dPredSum = currentElement |> Array.sumBy (fun x -> 
                                                                                    pq.[x].LeaveGroup idsParent
                                                                                    pq.[x].Top())
                                    let dCurrSum = currentElement |> Array.sumBy (fun x -> 
                                                                                    pq.[x].LeaveGroup currentGroupIDs
                                                                                    pq.[x].Top())
                                    {ParentDiss = dPredSum;
                                    CurrentDiss = dCurrSum}
                                else
                                    let dCurrSum = currentElement |> Array.sumBy (fun x -> 
                                                                                    pq.[x].LeaveGroup currentGroupIDs
                                                                                    pq.[x].Top())
                                    {element with CurrentDiss = dCurrSum})
                        else
                            group
                        )

            let newGainValue =
                gainValueFn nSet newGain

            let (downIter, best) =
                if (newGainValue < (fst bestSoFar)) then
                    (downIterLast+1, bestSoFar)
                else
                    (downIterLast, (newGainValue,currentConf))
        
            let cluster = currentConf.[clusterN]

            let outside =
                currentConf
                |> Array.concat
                |> Array.except cluster
            let (idOut, idIn, minDist) =
                outside
                |> Array.map (fun iOut -> Array.map (fun iIn -> (iOut, iIn, matrix.[iOut,iIn])) cluster)
                |> Array.concat
                |> Array.maxBy (fun (_, _, distance) -> distance) // find the closest item from outside,
            let (farItem, _, maxDist) =
                cluster
                |> Array.filter (fun iIn -> iIn<>idIn)
                |> Array.map (fun iIn -> (iIn, idIn, matrix.[iIn,idIn]))
                |> Array.minBy (fun (_, _, distance) -> distance) // find the farthest to the inside item within a cluster
            
            if (n >= 0) // size of the cluster as a limit for amount of iterations
                && (downIter<2) // no more than 2 down-iterations pro cluster
                && (maxDist) < (minDist) // there is a gain in swapping
                then
                    let idInsideItem = cluster |> Array.findIndex (fun id -> id = farItem)
                    let idOutsideCluster = currentConf |>  Array.findIndex (fun i -> i |> Array.exists (fun (ii) -> ii = idOut))
                    let idOutsideItem = Array.findIndex (fun id -> id = idOut) currentConf.[idOutsideCluster]
                    let temp = cluster.[idInsideItem]
                    currentConf.[clusterN].[idInsideItem] <- currentConf.[idOutsideCluster].[idOutsideItem]
                    currentConf.[idOutsideCluster].[idOutsideItem] <- temp
                    iterN <- iterN + 1
                    swap matrix clusterN currentConf (n-1) downIter best newGain [|clusterN;idOutsideCluster|] [|idInsideItem;idOutsideItem|]
                else best

        let fillSchemes (scheme: int []) : Map<string, Types.Item array> =
        
            let configurationIDs =
                scheme
                |> Array.mapFold (fun before n -> ([|before .. (before+n-1)|],(before+n))) 0
                |> fst
        
            let result =
                [|0 .. (scheme.Length-1)|]
                |> Array.fold (fun (gainBest, bestSoFar) clusterN ->
                    if (scheme.[clusterN]>1) then
                    
                        let gainComp =
                            configurationIDs
                            |> Array.map (fun idGroup ->
                                idGroup
                                |> Array.map (fun idElement ->
                                    let currentGroupIDs = idGroup |> Array.map (fun i -> initialGrouping.[i] |> snd) |> Array.concat  |> SSN.groupIDFn
                                    let currentElement = initialGrouping.[idElement] |> snd |> SSN.groupIDFn
                                    let dPredSum = currentElement |> Array.sumBy (fun x -> 
                                                                                    pq.[x].LeaveGroup idsParent
                                                                                    pq.[x].Top())
                                    let dCurrSum = currentElement |> Array.sumBy (fun x -> 
                                                                                    pq.[x].LeaveGroup currentGroupIDs
                                                                                    pq.[x].Top())
                                    {ParentDiss = dPredSum;
                                    CurrentDiss = dCurrSum}
                                )  )

                        swap matrixGainsForGroups clusterN bestSoFar scheme.[clusterN] 0 (0.,[||]) gainComp [||] [||]
                    else
                        (gainBest,bestSoFar)) (0., configurationIDs)
                |> snd
            
            result |> processInt
    
        let r =
            schemes
            |> Seq.rev
            |> List.ofSeq
            |> List.map fillSchemes
            |> Seq.ofList

        printfn "iterations done: %i" iterN

        r

    let kmeanSwapShufflePQ setN power f pq depth (map : Map<string,Types.Item array>) =
        let items = map |> Map.toArray
        [|for a in [1 .. power] do yield
                                        items
                                        //|> Array.shuffleFisherYates
                                        |> (fun x ->
                                                printfn "randomWalk %i" a
                                                kMeanSwapSaveGainPQ setN f pq depth x )
                                        |> Seq.toArray|]
        |> Array.stackVertical
        |> Array2D.toJaggedArray
        |> Array.map (fun a -> a |> Array.maxBy (fun m -> m |> Map.toArray |> Array.map (snd >> SSN.groupIDFn) |> gainstepcheckPQ f pq))


    
///////////////// old function: arrays, fixed initial values, no saving the best configuration

    let kMeanSwapOldsetK f matrixItems depth k (initialGrouping: (string*(Types.Item [])) []) : (Map<string,Types.Item []>) seq = 
    
        let parentGroup = initialGrouping |> Array.map (fun (k,p) -> p) |> Array.concat

        //let mutable iterN : int = 0

        let schemes =  SSN.schemeGenerator initialGrouping.Length k

        let gainDiff listA listB : float =
            let gFn (current: Types.Item []) = SSN.getStepGainFn f (SSN.groupIDFn current) (SSN.groupIDFn parentGroup) parentGroup.Length matrixItems
            (gFn (Array.append listA listB))

        let matrix =     
                let data =
                    initialGrouping 
                    |> Array.map snd
                let m = Array2D.zeroCreate (data.Length) (data.Length)
                for rowI in 0..data.Length-1 do
                    for colI in 0..rowI do
                        let tmp = if rowI=colI then nan else gainDiff data.[rowI] data.[colI]
                        m.[colI,rowI] <- tmp
                        m.[rowI,colI] <- tmp
                m

        let processInt (intListList: int [] []) =
            intListList
            |> Array.map (fun i ->
                                    if i.Length=1 then
                                        initialGrouping.[i.[0]]
                                    else
                                        Array.fold (fun (key,gr) ii ->
                                        let (binKey, value) = initialGrouping.[ii]
                                        ((sprintf "%s|%s" key binKey), Array.append value gr )
                                        ) ("",[||]) i)
            |> Array.map (fun (newBin,protA) ->
                (newBin, protA |> Array.map (fun prot ->
                        {prot with BinL = Array.append prot.BinL.[0..(depth)] [|newBin|]})))  
            |> Map.ofArray
    
        let rec swap (cluster: int []) (currentConf: int [] []) n = /// K-Mean-Swap for matrices
            
                let outside = 
                    currentConf
                    |> Array.concat
                    |> Array.except cluster
                let (idOut, idIn, minDist) =
                    outside
                    |> Array.map (fun i -> Array.map (fun j -> (i, j, matrix.[i,j])) cluster)
                    |> Array.concat
                    |> Array.maxBy (fun (outsideItem, insideItem, distance) -> distance) // find the closest item from outside, 
                let (farItem, closeItem, maxDist) =
                    cluster
                    |> Array.filter (fun iIn -> iIn<>idIn)
                    |> Array.map (fun i -> (i, idIn, matrix.[i,idIn])) 
                    |> Array.minBy (fun (farItem, closeItem, distance) -> distance) // find the farthest to the inside item within a cluster

                if (n >= 0)  then 
                    if (maxDist) < (minDist) then
                        let idInsideItem = cluster |> Array.findIndex (fun id -> id = farItem) 
                        let idOutsideCluster = currentConf |>  Array.findIndex (fun i -> i |> Array.exists (fun (ii) -> ii = idOut))
                        let idOutsideItem = Array.findIndex (fun id -> id = idOut) currentConf.[idOutsideCluster]
                        let temp = cluster.[idInsideItem]
                        cluster.[idInsideItem] <- currentConf.[idOutsideCluster].[idOutsideItem]
                        currentConf.[idOutsideCluster].[idOutsideItem] <- temp
                        //iterN <- iterN + 1
                        swap cluster currentConf (n-1) 
                //printfn "quit swap loop"

        let fillSchemes (scheme: int []) : Map<string, Types.Item []> =
    
            let mutable configurationIDs = 
                scheme
                |> Array.mapFold (fun before n -> ([|before .. (before+n-1)|],(before+n))) 0
                |> fst

            configurationIDs
            |> Array.iter (fun cluster -> if (cluster.Length>1) then swap cluster configurationIDs cluster.Length ) // not saving the best
        
            let res =
                configurationIDs
                |> processInt
            res

        //let r =                                   if want to see iteration count
        //    schemes
        //    |> Seq.rev
        //    |> List.ofSeq
        //    |> List.map fillSchemes
        //    |> Seq.ofList
        //printfn "iterations done: %i" iterN
        //r

        schemes
        |> Seq.rev
        |> Seq.map fillSchemes
    
    let kmeanSwapShuffleOldsetK power f matrixItems depth k (map : Map<string,Types.Item []>) =
        let items = map |> Map.toArray
        //kMeanSwapOld f matrixItems depth items
        //|> Seq.toArray
        [|for a in [1 .. power] do yield
                                        items
                                        |> Array.shuffleFisherYates 
                                        |> (fun x -> 
                                                printfn "randomWalk %i" a
                                                kMeanSwapOldsetK f matrixItems depth k x)
                                        |> Seq.toArray|]
        |> Array.stackVertical
        |> Array2D.toJaggedArray
        |> Array.map (fun a -> a |> Array.maxBy (fun m -> m |> Map.toArray |> Array.map (snd >> SSN.groupIDFn) |> gainstepcheck f matrixItems))



/////////////////////////       

/// call the main function, 
/// data is an experimental dataset for one functional group,
/// setN is usually the size of the whole dataset
let applySSN data setN = SSN.createTree (SSN.getStepGainNodeSetnR setN) (KMeanSwapFunctions.kmeanSwapShuffle setN 1) PreClusterFunctions.clusterHier (None) SSN.SSN data
let applySSNold data setN = SSN.createTree (SSN.getStepGainNodeSetnR setN) (KMeanSwapFunctions.kmeanSwapShuffleOld 1) PreClusterFunctions.clusterHier (None) SSN.SSN data
let applySSNnew data setN = SSN.createTree (SSN.getStepGainNodeSetnR setN) (KMeanSwapFunctions.kmeanSwapShuffleNew 1) PreClusterFunctions.clusterHier (None) SSN.SSN data
let readMM data setN = SSN.createTree (SSN.getStepGainNodeSetnR setN) (KMeanSwapFunctions.kmeanSwapShuffle setN 1) PreClusterFunctions.clusterHier (None) SSN.MM data

let readMMpq data setN = SSN.createTreePQ (SSN.getStepGainNodeSetnR setN) (KMeanSwapFunctions.kmeanSwapShuffle setN 1) PreClusterFunctions.clusterHier (None) SSN.MM data
let applySSNpq data setN = SSN.createTreePQ (SSN.getStepGainNodeSetnR setN) (KMeanSwapFunctions.kmeanSwapShuffle setN 1) PreClusterFunctions.clusterHier (None) SSN.SSN data

let applySSNpqKMS data setN = SSN.createTreePQinKMS (SSN.getStepGainNodeSetnR setN) (KMeanSwapFunctions.kmeanSwapShufflePQ setN 1) PreClusterFunctions.clusterHier (None) SSN.SSN data

let applySSN_Reduction k_pre data setN = SSN.createTreeReducing k_pre (SSN.getStepGainNodeSetnR setN) (KMeanSwapFunctions.kmeanSwapShuffle setN 1) PreClusterFunctions.kmeanGroups None SSN.SSN data
let applySSN_Reduction_Random k_pre data setN = SSN.createTreeReducing k_pre (SSN.getStepGainNodeSetnR setN) (KMeanSwapFunctions.kmeanSwapShuffle setN 1) PreClusterFunctions.kmeanGroupsRandom None SSN.SSN data

let applySSN_Reduction_Hier k_pre data setN = SSN.createTreeReducing k_pre (SSN.getStepGainNodeSetnR setN) (KMeanSwapFunctions.kmeanSwapShuffle setN 1) PreClusterFunctions.clusterHierGroups None SSN.SSN data

