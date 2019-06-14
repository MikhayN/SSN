module Functions

open Types
open PQ

open FSharp.Stats
open FSharpAux

open FSharp.Plotly

open System

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

    let getStepGainFn' fn itemsChild itemsParent (nRoot: int) matrix : float =
        let dSum itemsToSum itemsWhereToFind matrix =
        
            let findMaxDistIn idCurrent (idGroup: int array) (matrix: float [,]) = 
                idGroup
                |> Array.fold (fun maxSoFar i -> max maxSoFar matrix.[i,idCurrent]) 0.
        
            itemsToSum
            |> Array.fold (fun acc i -> acc + (findMaxDistIn i itemsWhereToFind matrix)) 0.

        let dPredSum = dSum itemsChild itemsParent matrix
        let dCurrSum = dSum itemsChild itemsChild matrix
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
        let schemeArray = Array.zeroCreate k 
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

    /// call hierarchical clustering for singletons
    let clusterHier (k: int) weight (children: Types.Item list) =
        let clusters nClusters =
            children
            |> List.map (fun protein -> protein.dataL)
            |> ML.Unsupervised.HierarchicalClustering.generate (weightedEuclidean weight) (ML.Unsupervised.HierarchicalClustering.Linker.centroidLwLinker)
            |> ML.Unsupervised.HierarchicalClustering.cutHClust nClusters
            |> List.map (List.map (fun i -> children.[ML.Unsupervised.HierarchicalClustering.getClusterId i]))
        clusters k
        |> List.map (fun list ->
            let binName = list |> List.fold (fun acc x -> sprintf "%s|p%i" acc x.ID) "hc" 
            (binName,list |> Array.ofList))
        |> Map.ofList
    
    // Pure combinatoric task

    ///Set partitioning
    let rec partitions l =
        let rec insertReplace x = function
            | []             -> []
            | (y :: ys) as l ->
                ((x::y)::ys)::(List.map (fun z -> y::z) (insertReplace x ys))
        seq {
            match l with   
            | []   ->
                yield []
            | h::tail ->
                for part in (partitions tail) do
                    yield [h]::part
                    yield! insertReplace h part
        }

    /// apply lazy partition to a broken groups
    let partGroup depth (groups: Map<string,Types.Item []>) : Map<string,Types.Item []> seq =
    
        let rename (list: (string*(Types.Item [])) list) =
            if list.Length>1 then
                let newKey =
                    list 
                    |> List.fold (fun state (k,v) -> (sprintf "%s|%s" state k)) (sprintf "mix")  
                let newListValue = 
                    list 
                    |> List.map (fun (k,v) -> v) 
                    |> Array.concat      
                    |> Array.map (fun protein -> {protein with BinL= Array.append protein.BinL.[0 .. depth] [|newKey|]})   
                (newKey, newListValue)
            else
                list.[0]
        groups
        |> Map.toList
        |> partitions
        |> Seq.map (List.map (fun variation -> rename variation) >> Map.ofList )

    //
    ////// Add KMean-Swap to find optimal conformation within a scheme
    //

    /// create tree function with two modes: MM - just read and show original MapMan ontology;
    /// SSN - process the MM tree into optimal SSN structure
    /// gainFn - gain formula, 
    /// kmeanswapFn - function for swapping,
    /// clusterFn - function for pre-clustering in case of more than 50 singletons as leaves
    let createTree gainFn (weight: seq<float> option) (mode: Types.Mode) (rootGroup: Types.Item array) = 
        
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
                |MM_raw -> // raw MapMan ontology without leaves breaking
                    let map = 
                        nodeMembers 
                        |> Array.filter (fun i -> i.BinL.Length > depth+1)
                        |> (fun i -> 
                                match i with
                                |[||] -> 
                                    Map.empty
                                |_ ->
                                    i
                                    |> Array.groupBy (fun i -> i.BinL.[depth+1]) 
                                    |> Map.ofArray)
                    if map=Map.empty then
                        Map.empty
                    else 
                        map
                        |> Map.map (fun key nodes -> 
                                            let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrix
                                            loop nodes (depth+1) dPredSum')
                |MM -> // MapMan with broken leaves
                    if nodeMembers.Length=1 then
                        Map.empty
                    else 
                        (breakGroup nodeMembers depth)
                        |> Map.map (fun key nodes -> 
                                        let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrix
                                        (loop nodes (depth+1) dPredSum'))

                |SSN (kmeanswapFn) -> // with simplification overall
                    if (nodeMembers.Length=1) 
                            || (nodeMembers.Length=0) 
                            || (String.contains "mix" (String.Concat nodeMembers.[0].BinL)) 
                            || (String.contains "|" (String.Concat nodeMembers.[0].BinL))
                            || (String.contains "c" (String.Concat nodeMembers.[0].BinL)) then 
                        Map.empty
                    else 
                        (breakGroup nodeMembers depth) // get a map with original grouping as (subbin, Item)
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

                |SSN_pre (kmeanswapFn, clusterFn) -> // with simplification overall
                    if (nodeMembers.Length=1) 
                            || (nodeMembers.Length=0) 
                            || (String.contains "mix" (String.Concat nodeMembers.[0].BinL)) 
                            || (String.contains "c" (String.Concat nodeMembers.[0].BinL)) then 
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

                |SSN_combi -> // without simplification, pure combinatorics
                    if (nodeMembers.Length=1) 
                            || (nodeMembers.Length=0) 
                            || (String.contains "mix" (String.Concat nodeMembers.[0].BinL)) 
                            || (String.contains "c" (String.Concat nodeMembers.[0].BinL)) 
                            || (String.contains "|" (String.Concat nodeMembers.[0].BinL)) then 
                        Map.empty
                    else 
                        (breakGroup nodeMembers depth)
                        |> partGroup depth
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

                |SST_walk (kmeanKKZ, walkFn) -> // with kMean as a start point for gain walking
                    if (nodeMembers.Length=1) 
                            || (nodeMembers.Length=0) 
                            || (String.contains "|" (String.Concat nodeMembers.[0].BinL)) then 
                        Map.empty
                    else 
                        (breakGroup nodeMembers depth)
                        |> (fun x -> 
                                                
                            let singles = 
                                x
                                |> Map.fold (fun state key nodes ->
                                        let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrix
                                        state 
                                        |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)

                            walkFn kmeanKKZ depth matrix singles gainFn x 
                            |> Map.fold (fun state key nodes ->
                                match (singles.TryFind key) with
                                | None ->
                                    let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrix
                                    state 
                                    |> Map.add key (loop nodes (depth+1) dPredSum')
                                | Some x ->
                                    state 
                                    |> Map.add key x                                                                   
                            ) (Map.empty) )

            let confGain = confGainFn children
            {
            Member = nodeMembers;
            Children = 
                match mode with
                |MM -> 
                    children
                |MM_raw -> 
                    children
                |_ -> 
                    if  (confGain > stepGain) then 
                        children;
                    else 
                        Map.empty
            StepGain = stepGain; 
            ConfGain = (confGain, children |> Map.toList |> List.map fst);
            GroupGain = max stepGain confGain;
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

module KMeanSwapFunctions =
    
    let gainstepcheck f (matrix: float [,]) (conf: int array array)  =
        let parentGroup = conf |> Array.concat
        let gFn (current: int array) = SSN.getStepGainFn f current parentGroup parentGroup.Length matrix
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
                                        ) ("mix",[||]) i)
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
        kMeanSwapOld f matrixItems depth items
        |> Seq.toArray
        //[|for a in [1 .. power] do yield
        //                                items
        //                                //|> Array.shuffleFisherYates 
        //                                |> (fun x -> 
        //                                        printfn "randomWalk %i" a
        //                                        kMeanSwapOld f matrixItems depth (x))
        //                                |> Seq.toArray|]
        //|> Array.stackVertical
        //|> Array2D.toJaggedArray
        //|> Array.map (fun a -> a |> Array.maxBy (fun m -> m |> Map.toArray |> Array.map (snd >> SSN.groupIDFn) |> gainstepcheck f matrixItems))

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
                                        ) ("mix",[||]) i)
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
                                        |> Array.shuffleFisherYates 
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

    let kMeanSwapSaveGain nSet f matrixSingletons depth (initialGrouping: (string*(Types.Item array)) []) : (Map<string,Types.Item array>) seq =
    
        let parentGroup = initialGrouping |> Array.map (fun (k,p) -> p) |> Array.concat

        let mutable iterN : int = 0

        let schemes =  seq [2..initialGrouping.Length] |> Seq.map (SSN.schemeGenerator initialGrouping.Length) |> Seq.concat

        let gainDiff listA listB : float =
            let gFn (current: Types.Item array) = SSN.getStepGainFn f (SSN.groupIDFn current) (SSN.groupIDFn parentGroup) parentGroup.Length matrixSingletons
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
                                        ) ("mix",[||]) i)
            |> Array.map (fun (newBin,protA) ->
                (newBin, protA |> Array.map (fun prot ->
                        {prot with BinL = Array.append prot.BinL.[0..(depth)] [|newBin|]}))) // updating the bin labels            
            |> Map.ofArray
    
        let gainValueFn (gainA: GainComponents [] []) =
            gainA
            |> Array.fold (fun acc group ->
                let nC = float group.Length
                let nR = float nSet
                let deltaSpec = -((nC/nR)*log2(nC/nR)+((nR-nC)/nR)*log2((nR-nC)/nR))
                acc + (group |> Array.sumBy (fun x -> x.ParentDiss-x.CurrentDiss))*deltaSpec
                ) 0.

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
                gainValueFn newGain

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

        printfn "iterations done: %i" iterN

        r

    let kmeanSwapShuffle setN power f matrixItems depth (map : Map<string,Types.Item array>) =
        let items = map |> Map.toArray
        [|for a in [1 .. power] do yield
                                        items
                                        //|> Array.shuffleFisherYates
                                        |> (fun x ->
                                                printfn "randomWalk %i" a
                                                kMeanSwapSaveGain setN f matrixItems depth x )
                                        |> Seq.toArray|]
        |> Array.stackVertical
        |> Array2D.toJaggedArray
        |> Array.map (fun a -> a |> Array.maxBy (fun m -> m |> Map.toArray |> Array.map (snd >> SSN.groupIDFn) |> gainstepcheck f matrixItems))

module Clustering =
    
    
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

    /// Kmean for groups with KKZ centroid init
    let kmeanGroupsKKZ (k: int) (children: Map<string,Types.Item []> ) =

        let data = children |> Map.toArray |> Array.map (fun (s,ar) -> (s,ar |> Array.map (fun p -> p.dataL) |> MatrixTopLevelOperators.matrix))

        let clusters = ML.Unsupervised.IterativeClustering.compute distMatrix (initCgroupsKKZ) updateCentroid data k

        data
        |> Array.map (fun list -> (clusters.Classifier list |> fst),list )
        |> Array.groupBy (fst)
        |> Array.map (fun (_,list) -> 
            list 
            |> Array.map (fun (_,(bin,_)) -> (bin, children |> Map.find bin))
            )

module Walk =

    let reflectStateID (matrixA) =
        matrixA 
        |> List.distinct
        |> List.map (fun i -> 
            i 
            |> List.indexed 
            |> List.filter (fun (_,x) -> x=1 ) 
            |> List.map (fun (i,_) -> i)
            |> List.toArray
            )
        |> List.toArray
    
    let gainFnn (data: Item [] []) (singleGG: float []) matrixSingles fn (x: int [] []) = 
            x
            |> Array.sumBy (fun ids ->
                    if ids.Length=1 then
                        singleGG.[ids.[0]]
                    else    
                        let itemsChild = ids |> Array.map (fun i -> data.[i]) |> Array.concat |> SSN.groupIDFn
                        let itemsParent = data |> Array.concat |> SSN.groupIDFn 
                        SSN.getStepGainFn fn itemsChild itemsParent itemsParent.Length matrixSingles
            )

    let matrixG_from_matrixA (data: Item [] []) (singleGG: float []) matrixSingles fn matrixA =
        let clusterMA = matrixA |> List.distinct
        let pairArray = Array.allPairs [|0 .. (data.Length-1)|] [|0 .. (clusterMA.Length)|] 
        let m = JaggedArray.zeroCreate data.Length (clusterMA.Length+1)
        pairArray
        |> Array.iter (fun (i,J) -> // that can be substituted by for i=0 to data.Length do for J=0 to (clusterMA.Length+1) do
            
            let matrixAA = matrixA |> List.map (List.toArray) |> List.toArray

            let A = (matrixAA.[i] |> Array.indexed |> Array.filter (fun (i,x) -> x=1) |> Array.map fst |> Array.toList) // a was in A
            
            let B = 
                if J<clusterMA.Length then
                    (clusterMA.[J] |> List.indexed |> List.filter (fun (i,x) -> x=1) |> List.map fst) // b was in B
                else 
                    if A.Length=1 then
                        []
                    else
                        [i]

            if (A=B) || (B.IsEmpty) then 
                m.[i].[J] <- (0., [])
            
            else
                for ii in A do   // move 'a' out of 'A', but leave it stay with itself
                        matrixAA.[i].[ii] <- 0
                        matrixAA.[ii].[i] <- 0
                matrixAA.[i].[i] <- 1

                for jj in B do   // move a in B
                        matrixAA.[i].[jj] <- 1
                        matrixAA.[jj].[i] <- 1
                let newMA = matrixAA |> Array.map (Array.toList) |> Array.toList
                        
                let stateIDs = reflectStateID newMA
                let gain =
                    stateIDs |> gainFnn data singleGG matrixSingles fn

                m.[i].[J] <- (gain, newMA)
            )
        m


    //type QDictionary_GValue =
    //    {
    //    NextStepGain: float
    //    NextStepState: int list list
    //    mutable MaxGain: float
    //    } 

    let walkingFn dN sN kmeanKKZ depth matrixSingletons (singles: Map<string,Node<string,Item>>) gainFn (data: Map<string, Item []>) = 
    
        let mutable qDictionary: Map<(int list list),((float*(int list list)) [] [])> = Map.empty

        let dataGroupsA = data |> Map.toArray

        let singleGG =
            dataGroupsA
            |> Array.map (fun (bin,_) -> 
                                let node = singles |> Map.find bin 
                                node.GroupGain)

        /// input data: groups of items
        let dataGroups = dataGroupsA |> Array.map snd

        let rename (list: (string*(Types.Item [])) []) =
                if list.Length>1 then
                    let newKey =
                        list 
                        |> Array.fold (fun state (k,v) -> (sprintf "%s|%s" state k)) (sprintf "mix")  
                    let newListValue = 
                        list 
                        |> Array.map (fun (k,v) -> v) 
                        |> Array.concat      
                        |> Array.map (fun protein -> {protein with BinL= Array.append protein.BinL.[0 .. depth] [|newKey|]})   
                    (newKey, newListValue)
                else
                    list.[0]    

        let superFunctionTestG (data: Item [] []) (singleGG: float []) fn matrixSingles initConf  =

            let initialConfig = 
                initConf
                |> Array.map 
                    (Array.map (fun (label,_) -> 
                        dataGroupsA 
                        |> Array.findIndex (fun (l,_) -> l=label) )
                    )

            //let fileLogInitState = sprintf "Initial state: %A"  (initConf |> Array.map (Array.map fst))
            //File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogInitState])

            /// adjustency matrix (initialize with initial configuration) 
            let matrixA =
                let m = JaggedArray.zeroCreate data.Length data.Length
                for i=0 to (data.Length-1) do
                    let cluster = initialConfig |> Array.find (fun x -> x |> Array.contains i)
                    for j=0 to (data.Length-1) do
                        if (cluster |> Array.contains j) then   
                            m.[i].[j] <- 1
                        else
                            m.[i].[j] <- 0
                m
                |> Array.map (Array.toList)
                |> Array.toList

            //let mutable qDictionary': Map<(int list list),(QDictionary_GValue [] [])> = Map.empty
            //((qDictionary'.Item matrixA).[0].[0]).MaxGain <- 0. 

            let initStateIDs = reflectStateID matrixA
            let initialState = (gainFnn data singleGG matrixSingles fn initStateIDs, initStateIDs)

            //let fileLogInit = sprintf "0\t0\tnan\tnan\t0.\t%f" (fst initialState)
            //File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogInit])

            //let fileLog = sprintf "0\t0\t0.\t%f" (fst initialState)
            //File.AppendAllLines((sprintf "%s%s.txt" pathLog fileLogName), [fileLog])

            let clusterMA = matrixA |> List.distinct
            let pairArray' = Array.allPairs [|0 .. (data.Length-1)|] [|0 .. (clusterMA.Length)|] 

            let matrixG_origin = matrixG_from_matrixA data singleGG matrixSingles fn matrixA
        
            qDictionary <- (qDictionary.Add (matrixA, matrixG_origin))

            let pq_origin =
                let n = pairArray'.Length
                let p = MaxIndexPriorityQueue<float>(n)
                for id=0 to n-1 do 
                    let (i,ii) = pairArray'.[id]
                    if (fst matrixG_origin.[i].[ii])>0. then p.Insert id (fst matrixG_origin.[i].[ii]) // load all calculated gains
                p

            let rec loop iStep (mA: int list list) (pairArray: (int*int) []) (mG: (float*(int list list)) [] []) (pq: MaxIndexPriorityQueue<float>) (moved: int []) =
        
                let gainCurrent = gainFnn data singleGG matrixSingles fn (reflectStateID mA)
        
                let mutable countDirections = 0
            
                seq [ while 
                    (pq.Length>0) 
                    && (pq.Top()>0.) 
                    && (countDirections<dN) && (iStep<sN) do // (pq.Top() > gainCurrent) && (countDirections<3) && (iStep<6)
                
                        countDirections <- countDirections + 1
                    
                        // order represents the moving: a - moved element, b - target cluster
                                     

                        //let fileLogStep = sprintf "%i\t%i\t%i\t%i\t%f\t%f" iStep countDirections a b gainCurrent (pq.Top())
                        //File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogStep])

                        let (a,b) = 
                            while (pq.Length>0) 
                                && (qDictionary |> Map.containsKey (snd mG.[fst pairArray.[pq.TopIndex()]].[snd pairArray.[pq.TopIndex()]])) do
                                    pq.Pop() |> ignore

                            if pq.Length=0 then 
                                (-1,-1)
                            else
                                pairArray.[pq.TopIndex()]
                        
                        if (a,b)=(-1,-1) then // if the pq is empty and no new states are found
                            yield! []
                        else

                            let mA_new = snd mG.[a].[b]

                            //let fileStep = sprintf "%i\t%i\t%f\t%f" iStep countDirections gainCurrent (pq.Top())
                            //File.AppendAllLines((sprintf "%s%s.txt" pathLog fileLogName), [fileStep])

                            pq.Pop() |> ignore // pq will be used for other directiones in while loop

                            // find all values in mG with the same state and exclude possibility to go there again 
                            // by adding all a's in moved (don't change mG!) and removing duplicate states from pq
                            let all_a = 
                                mG 
                                |> Array.indexed 
                                |> Array.filter (fun (_, vL) -> (vL |> Array.contains mG.[a].[b]) )
                                |> Array.map (fun (i,vl) ->
                                            let jjj = vl |> Array.findIndex (fun v -> v = mG.[a].[b])
                                            pq.TryRemove ((i*mG.[0].Length)+jjj)
                                            i )

                        //if (qDictionary |> Map.containsKey mA_new) then
                        
                        //    //let fileLogState = sprintf "%A was already visited, no step further" (reflectStateID mA_new)
                        //    //File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogState]) 
                        
                        //    yield! [] // how to get rid of the unnecessary empty lists? 
                        //else
                    
                            //let fileLogState = sprintf "%A" (reflectStateID mA_new)
                            //File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogState]) 

                            let clusterMA_new = mA_new |> List.distinct
                            let pairArrayNew = Array.allPairs [|0 .. (data.Length-1)|] [|0 .. (clusterMA_new.Length)|] 

                            let matrixG = matrixG_from_matrixA data singleGG matrixSingles fn mA_new
                        
                            qDictionary <- (qDictionary.Add (mA_new, matrixG))
                        
                            let pq_new = 
                                let n = pairArrayNew.Length
                                let p = MaxIndexPriorityQueue<float>(n)
                                for j=0 to n-1 do 
                                    let (i,ii) = pairArrayNew.[j]
                                    let gain = fst matrixG.[i].[ii]
                                    if gain>0. then
                                        p.Insert j (gain) // load all gains except of redundant
                                p
                        
                            let new_moved = Array.append all_a moved |> Array.distinct

                            new_moved
                            |> Array.iter (fun i ->
                                let indices = [|(i * matrixG.[0].Length) .. (i * (matrixG.[0].Length) + matrixG.[0].Length - 1)|]
                                pq_new.TryRemoveGroup indices
                            )

                            let configuration = reflectStateID mA_new
                            let gain = gainFnn data singleGG matrixSingles fn (configuration)
                            let stats = (gain, configuration) 
                            
                            yield (stats)
                            yield! loop (iStep+1) mA_new pairArrayNew matrixG pq_new new_moved
                ]
    
            Seq.appendSingleton (loop 1 matrixA pairArray' matrixG_origin pq_origin [||]) initialState 

        (seq [singleGG |> Array.sum, Array.init data.Count (fun i -> [|i|])]) :: 
            ([2 .. (data.Count-1)] 
            |> List.map (fun i -> 

                //let fileLogKMEAN = sprintf "Kmean with k=%i"  i
                //File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogKMEAN])

                data
                |> kmeanKKZ i
                |> superFunctionTestG dataGroups singleGG gainFn matrixSingletons)
            )
        |> Seq.ofList
        |> Seq.concat           
        |> Seq.maxBy (fst)
        |> snd
        |> Array.map (fun groupIDs ->     
            groupIDs 
            |> Array.map (fun groupID -> 
                dataGroupsA.[groupID])
            |> rename )
        |> Map.ofArray

/////////////////////////       

/// call the main function, 
/// data is an experimental dataset for one functional group,
/// setN is usually the size of the whole dataset
let readMM_raw setN data  = SSN.createTree (SSN.getStepGainNodeSetnR setN)  (None) Types.Mode.MM_raw data
let readMM setN data  = SSN.createTree (SSN.getStepGainNodeSetnR setN)  (None) Types.Mode.MM data
let applySSN setN data = SSN.createTree (SSN.getStepGainNodeSetnR setN)  (None) (Types.Mode.SSN_pre ((KMeanSwapFunctions.kmeanSwapShuffle setN 1), SSN.clusterHier)) data
let applySSNcombi setN data = SSN.createTree (SSN.getStepGainNodeSetnR setN)  (None) Types.Mode.SSN_combi data
let applySSNold setN data = SSN.createTree (SSN.getStepGainNodeSetnR setN) (None) (SSN_pre ((KMeanSwapFunctions.kmeanSwapShuffleOld 1), SSN.clusterHier)) data
let applySST_walk setN data = SSN.createTree (SSN.getStepGainNodeSetnR setN) None (SST_walk (Clustering.kmeanGroupsKKZ, (Walk.walkingFn 1 5))) data

let asyncApplySSN setN data = async {return (applySST_walk setN data)}





