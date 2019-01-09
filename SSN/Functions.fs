module Functions

open Types

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
                            || (String.contains "c" (String.Concat nodeMembers.[0].BinL)) then 
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

            let confGain = confGainFn children
            {
            Member = nodeMembers;
            Children = 
                match mode with
                |MM -> 
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

/////////////////////////       

/// call the main function, 
/// data is an experimental dataset for one functional group,
/// setN is usually the size of the whole dataset
let readMM setN data  = SSN.createTree (SSN.getStepGainNodeSetnR setN)  (None) Types.Mode.MM data
let applySSN setN data = SSN.createTree (SSN.getStepGainNodeSetnR setN)  (None) (Types.Mode.SSN_pre ((KMeanSwapFunctions.kmeanSwapShuffle setN 1), SSN.clusterHier)) data
let applySSNcombi setN data = SSN.createTree (SSN.getStepGainNodeSetnR setN)  (None) Types.Mode.SSN_combi data
let applySSNold setN data = SSN.createTree (SSN.getStepGainNodeSetnR setN) (None) (SSN_pre ((KMeanSwapFunctions.kmeanSwapShuffleOld 1), SSN.clusterHier)) data


