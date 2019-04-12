﻿//#r @"c:\Users\mikha\source\repos\mathnet-numerics\src\Numerics\bin\Debug\netstandard2.0\MathNet.Numerics.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\netstandard.dll"

//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\BioFSharp.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\BioFSharp.IO.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharpAux.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharpAux.IO.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\Newtonsoft.Json.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharp.Stats.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharp.Plotly.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharpGephiStreamer.dll"

//#load "Types.fs"
//#load "PQ.fs"
//#load "Functions.fs"
//#load "FunctionsExp.fs"
//#load "GePhi.fs"
//#load "TestData.fs"
//#load "SOM.fs"
//#load "Auxilliary.fs"
//#load "Plots.fs"

//open System 
//open System.IO 
//open FSharpAux
//open FSharp.Plotly
//open FSharp.Stats

//open Functions.SSN
//open TestData
//open GePhi
//open Types
//open PQ
//open Auxilliary
//open Plots
//open FunctionsExp

//#time

///// file subname for the trial # change each time for the new file!!!!! ##############################################################################
//let fileSubName = "ProteinDataSet_Path6"

///// file for tracking a path during the whole treeCreate. 
//let header = "StepN\tDirectionN\tmovedID\tMovedTo\tGainBefore\tGainNow" 
//let pathLOG = @"c:\Users\mikha\Work-CSB\Projects\SSN\results\walkResults\PathLog_"
//File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [header])


//let mutable stepCount = 0

//let pairArrayFn (groupL: int) = 
//    [|0 .. (groupL-1)|] 
//    |> Array.allPairs [|0 .. (groupL-1)|] 
//    |> Array.filter (fun (i,j) -> i<>j)

//let reflectStateID (matrixA) =
//    matrixA 
//    |> List.distinct
//    |> List.map (fun i -> 
//        i 
//        |> List.indexed 
//        |> List.filter (fun (_,x) -> x=1 ) 
//        |> List.map (fun (i,_) -> i)
//        |> List.toArray
//        )
//    |> List.toArray
    
////let reflectStateG (matrixA) (data: Item [] []) =
////    matrixA 
////    |> List.distinct
////    |> List.map (fun i -> 
////        i 
////        |> List.indexed 
////        |> List.filter (fun (_,x) -> x=1 ) 
////        |> List.map (fun (i,_) -> data.[i])
////        |> List.toArray
////        |> Array.concat
////        )
////    |> List.toArray
    

///// loop for path walking with nDirections iterations and nSteps in depth, applying f on each step
//let superFunctionTestG (data: Item [] []) (singleGG: float []) fn matrixSingles matrixA  =
    
//    let mutable qDictionary: Map<(int list list),((float*(int list list)) [] [])> = Map.empty

//    let initStateIDs = reflectStateID matrixA
    
//    let gainFnn (x: int [] []) = 
//        x
//        |> Array.sumBy (fun ids ->
//                if ids.Length=1 then
//                    singleGG.[ids.[0]]
//                else    
//                    getStepGainFn fn (ids |> Array.map (fun i -> data.[i]) |> Array.concat |> groupIDFn) (data |> Array.concat |> groupIDFn ) (data |> Array.concat |> Array.length) matrixSingles
//        )

//    let initialState = (gainFnn initStateIDs, initStateIDs)

//    let fileLogInit = sprintf "0.\t0.\tnan\tnan\t0.\t%f" (fst initialState)
//    File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogInit])

//    let clusterMA = matrixA |> List.distinct
//    let pairArray = Array.allPairs [|0 .. (data.Length-1)|] [|0 .. (clusterMA.Length)|] 

//    let matrixG_origin =
//        let m = JaggedArray.zeroCreate data.Length (clusterMA.Length+1)
//        pairArray
//        |> Array.iter (fun (i,J) ->     
            
//            let matrixAA = matrixA |> List.map (List.toArray) |> List.toArray

//            let A = (matrixAA.[i] |> Array.indexed |> Array.filter (fun (i,x) -> x=1) |> Array.map fst |> Array.toList) // a was in A
            
//            let B = 
//                if J<clusterMA.Length then
//                    (clusterMA.[J] |> List.indexed |> List.filter (fun (i,x) -> x=1) |> List.map fst) // b was in B
//                else 
//                    if A.Length=1 then
//                        []
//                    else
//                        [i]

//            if (A=B) || (B.IsEmpty) then 
//                m.[i].[J] <- (0., [])
            
//            else
//                for ii in A do   // move 'a' out of 'A', but leave it stay with itself
//                        matrixAA.[i].[ii] <- 0
//                        matrixAA.[ii].[i] <- 0
//                matrixAA.[i].[i] <- 1

//                for jj in B do   // move a in B
//                        matrixAA.[i].[jj] <- 1
//                        matrixAA.[jj].[i] <- 1
//                let newMA = matrixAA |> Array.map (Array.toList) |> Array.toList
                        
//                let stateIDs = reflectStateID newMA
//                let gain =
//                    stateIDs |> gainFnn

//                m.[i].[J] <- (gain, newMA)
//            )
//        m

//    qDictionary <- (qDictionary.Add (matrixA, matrixG_origin))

//    let pq_origin =
//        let n = pairArray.Length
//        let p = MaxIndexPriorityQueue<float>(n)
//        for id=0 to n-1 do 
//            let (i,ii) = pairArray.[id]
//            if (fst matrixG_origin.[i].[ii])>0. then p.Insert id (fst matrixG_origin.[i].[ii]) // load all calculated gains
//        p

//    let rec loop iStep (mA: int list list) (pairArray: (int*int) []) (mG: (float*(int list list)) [] []) (pq: MaxIndexPriorityQueue<float>) (moved: int list) =
        
//        let gainCurrent = gainFnn (reflectStateID mA)
        
//        let mutable countDirections = 0
//        //printfn "compare pq.Top:%f vs gain:%f" (pq.Top()) gainCurrent
//        seq [ while (pq.Top()>0.) && (pq.Top() > gainCurrent) && (countDirections<3) && (iStep<15) do
                
//                countDirections <- countDirections + 1

//                let pairID = pq.TopIndex()

//                // order represents the moving: a - moved element, b - target cluster
//                let (a,b) = pairArray.[pairID] 

//                let fileLogStep = sprintf "%i\t%i\t%i\t%i\t%f\t%f" iStep countDirections a b gainCurrent (pq.Top())
//                File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogStep])

//                pq.Pop() |> ignore           // pq will be used for other directiones in while loop

//                let mA_new = snd mG.[a].[b]
                    
//                let fileLogState = sprintf "%A" (reflectStateID mA_new)
//                File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogState]) 

//                let clusterMA_new = mA_new |> List.distinct
//                let pairArrayNew = Array.allPairs [|0 .. (data.Length-1)|] [|0 .. (clusterMA_new.Length)|] 

//                let matrixG =
//                    let mTemp = 
//                        if (qDictionary |> Map.containsKey mA_new) then
//                            qDictionary |> Map.find mA_new
//                        else   
//                            let m = JaggedArray.zeroCreate data.Length (clusterMA_new.Length+1)
//                            pairArrayNew
//                            |> Array.iter (fun (i,J) ->     
//                                let matrixAA = mA_new |> List.map (List.toArray) |> List.toArray

//                                let A = (matrixAA.[i] |> Array.indexed |> Array.filter (fun (i,x) -> x=1) |> Array.map fst |> Array.toList) // a was in A
            
//                                let B = 
//                                    if J<clusterMA_new.Length then
//                                        (clusterMA_new.[J] |> List.indexed |> List.filter (fun (i,x) -> x=1) |> List.map fst) // b was in B
//                                    else 
//                                        if A.Length=1 then
//                                            []
//                                        else
//                                            [i]

//                                if (A=B) || (B.IsEmpty) then 
//                                    m.[i].[J] <- (0., [])
     
//                                else
//                                    for ii in A do   // move 'a' out of 'A', but leave it stay with itself
//                                            matrixAA.[i].[ii] <- 0
//                                            matrixAA.[ii].[i] <- 0
//                                    matrixAA.[i].[i] <- 1

//                                    for jj in B do   // move a in B
//                                            matrixAA.[i].[jj] <- 1
//                                            matrixAA.[jj].[i] <- 1
//                                    let newMA = matrixAA |> Array.map (Array.toList) |> Array.toList
                        
//                                    let stateIDs_new = reflectStateID newMA
//                                    let gain_new =
//                                        stateIDs_new |> gainFnn
//                                    m.[i].[J] <- (gain_new, newMA)
//                                )
                            
//                            m
//                    qDictionary <- (qDictionary.Add (mA_new, mTemp))
//                    (a::moved) |> List.iter (fun xa -> mTemp.[xa] <- Array.create (clusterMA_new.Length+1) (0., []))
//                    mTemp

//                let pq_new = 
//                    let n = pairArrayNew.Length
//                    let p = MaxIndexPriorityQueue<float>(n)
//                    for j=0 to n-1 do 
//                        let (i,ii) = pairArrayNew.[j]
//                        let gain = fst matrixG.[i].[ii]
//                        if gain>0. then
//                            p.Insert j (gain) // load all gains except of (a -> a)
//                    p

//                let configuration = reflectStateID mA_new
//                let gain = gainFnn (configuration)
//                let stats = (gain, configuration) //, iDirection, iStep)

//                stepCount <- stepCount + 1

//                if (pq_new.Top() <= gain) then
//                    yield (stats)
//                else
//                    yield (stats)
//                    yield! loop (iStep+1) mA_new pairArrayNew matrixG pq_new (a::moved)

//        ]
    
//    let r = loop 1 matrixA pairArray matrixG_origin pq_origin []
//    Seq.appendSingleton r initialState
//     //|> Seq.toList

//// adjust it to following without explicit nSteps and nDir
//let plotGainWalk (list: ('a * 'b * int * int) list ) =
//    let color = function
//        |1 -> colorBlue
//        |2 -> colorBrightBlue
//        |3 -> colorGreen
//        |4 -> colorOrange
//        |5 -> colorYellow
//        |_ -> colorGray
    
//    let edges = 
//        list 
//        |> List.tail
//        |> List.mapi ( fun i (_,_,_,dCurrent) -> 
                                
//                                let source = 
//                                    (list 
//                                    |> List.cutAfterN (i+1)
//                                    |> fst
//                                    |> List.rev
//                                    |> List.tryFindIndex (fun (_,_,_,dBefore) -> dBefore < dCurrent)) 
//                                let target = i+1
//                                match source with
//                                |None -> (0,target)
//                                |Some x -> (i - x,target)
//                                )
//    let lines = 
//        edges
//        |> List.map (fun (i,j) -> 
//                        let (gain1,_,_,step1) = list.[i]
//                        let (gain2,_,_,step2) = list.[j]
//                        Chart.Line ([(step1,gain1);(step2,gain2)], Color = colorGray)
//                    )
//    let points =
//        list
//        |> List.map (fun (gain,_,direction,step) -> 
//                        Chart.Point ([(step,gain)], Labels = [sprintf "%i" direction] ) |> Chart.withMarkerStyle (Size=10, Color=color direction)
//                        )

//    [lines; points] |> List.concat |> Chart.Combine |> Chart.Show

////plotGainWalk dataToPlot



//// incorporate the pathWalking into the createTree

//type ClusterFn = int -> Map<string, Item []> -> (string * (Item [])) [] []

//type WalkingFn = int -> float [,] -> Map<string,Node<string,Item>> -> (float -> float -> int -> int -> float) -> Map<string, Item []> -> (string * (Item [])) [] [] -> Map<string, Item []>

//type Mode =
//    |MM                                     // MapMan with broken leaves
//    |SSN_combi                              // without simplification, pure combinatorics
//    |SSN_walk of (ClusterFn*WalkingFn)      // with kMean as a start point for gain walking


//let createTreeWalking gainFn (weight: seq<float> option) (mode: Mode) (rootGroup: Types.Item array) = 
        
//    let nRoot = rootGroup.Length

//    let matrixSingletons = 
//        rootGroup
//        |> distMatrixWeightedOf SSN.distanceMatrixWeighted weight

//    // calculation for one node    
//    let rec loop (nodeMembers: Types.Item array) depth dPredSum =
        
//        /// sum of max dist within a node 
//        let dCurrSum = dSumFn (groupIDFn nodeMembers) (groupIDFn nodeMembers) matrixSingletons

//        /// to calc step from parent to current node
//        let stepGain = (gainFn dCurrSum dPredSum nodeMembers.Length nRoot)

//        let children = 
//            match mode with
//            |MM -> // MapMan with broken leaves
//                if nodeMembers.Length=1 then
//                    Map.empty
//                else 
//                    (breakGroup nodeMembers depth)
//                    |> Map.map (fun key nodes -> 
//                                    let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrixSingletons
//                                    (loop nodes (depth+1) dPredSum'))


//            |SSN_combi -> // without simplification, pure combinatorics
//                if (nodeMembers.Length=1) 
//                        || (nodeMembers.Length=0) 
//                        || (String.contains "|" (String.Concat nodeMembers.[0].BinL)) then 
//                    Map.empty
//                else 
//                    (breakGroup nodeMembers depth)
//                    |> partGroup depth
//                    |> Seq.fold (fun (singles,best) i -> 
//                        let newNodes = 
//                            if singles=Map.empty then
//                                i
//                                |> Map.fold (fun state key nodes ->
//                                    let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrixSingletons
//                                    state 
//                                    |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)
//                            else
//                                i
//                                |> Map.fold (fun state key nodes ->
//                                    match (singles.TryFind key) with
//                                    | None ->
//                                        let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrixSingletons
//                                        state 
//                                        |> Map.add key (loop nodes (depth+1) dPredSum')
//                                    | Some x ->
//                                        state 
//                                        |> Map.add key x                                                                   
//                                ) (Map.empty)
                                                                    
//                        let best' =
//                            if (confGainFn newNodes) > (confGainFn best) then  // compare configuration gains to get the best
//                                newNodes  
//                            else 
//                                best
//                        if (singles = Map.empty) then
//                            (newNodes, best')
//                        else
//                            (singles, best')
                                            
//                    ) (Map.empty,Map.empty) // here as state should be this singles (first) saved and optimal conf
//                    |> snd

//            |SSN_walk (kmeanKKZ, walkFn) -> // with kMean as a start point for gain walking
//                if (nodeMembers.Length=1) 
//                        || (nodeMembers.Length=0) 
//                        || (String.contains "|" (String.Concat nodeMembers.[0].BinL)) then 
//                    Map.empty
//                else 
//                    (breakGroup nodeMembers depth)
//                    |> (fun x -> 
//                        let nMax = x.Count
                        
//                        let fileLogData = sprintf "for elements:%A" (x |> Map.toArray |> Array.map fst)
//                        File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogData])

//                        let singles = 
//                            x
//                            |> Map.fold (fun state key nodes ->
//                                    let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrixSingletons
//                                    state 
//                                    |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)

//                        //printfn "singleGG %A" singles

                        
//                        [2 .. (nMax-1)] 
//                        |> List.rev
//                        |> List.map (fun i -> 

//                            let fileLogKMEAN = sprintf "Kmean with k=%i"  i
//                            File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogKMEAN])

//                            x
//                            |> kmeanKKZ i
//                            |> walkFn depth matrixSingletons singles gainFn x) // here the first in seq should be x, then everything else
                        
//                        |> Seq.fold (fun (singles,best) i -> 
//                            let newNodes = 
//                                if singles=Map.empty then
//                                    i
//                                    |> Map.fold (fun state key nodes ->
//                                        let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrixSingletons
//                                        state 
//                                        |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)
//                                else
//                                    i
//                                    |> Map.fold (fun state key nodes ->
//                                        match (singles.TryFind key) with
//                                        | None ->
//                                            let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrixSingletons
//                                            state 
//                                            |> Map.add key (loop nodes (depth+1) dPredSum')
//                                        | Some x ->
//                                            state 
//                                            |> Map.add key x                                                                   
//                                    ) (Map.empty)
                                                                    
//                            let best' =
//                                if (confGainFn newNodes) > (confGainFn best) then  // compare configuration gains to get the best
//                                    newNodes  
//                                else 
//                                    best
//                            (singles, best')
                                            
//                        ) (singles,singles) )// here as state should be this singles (first) saved and optimal conf
//                    |> snd

//        let confGain = confGainFn children
//        {
//        Member = nodeMembers;
//        Children = 
//            match mode with
//            |MM -> 
//                children
//            |_ -> 
//                if  (confGain > stepGain) then 
//                    children;
//                else 
//                    Map.empty
//        StepGain = stepGain; 
//        ConfGain = (confGain, children |> Map.toList |> List.map fst);
//        GroupGain = max stepGain confGain;
//        }
    
//    loop rootGroup 0 0.

//let walkingFn depth matrixSingletons (singles: Map<string,Node<string,Item>>) gainFn (data: Map<string, Item []>) (initConfig: (string * (Item [])) [] [])  = 
    
//    let fileLogInitState = sprintf "Initial state: %A"  (initConfig |> Array.map (Array.map fst))
//    File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogInitState])

//    let dataGroupsA = data |> Map.toArray

//    let singleGG =
//        dataGroupsA
//        |> Array.map (fun (bin,_) -> 
//                            let node = singles |> Map.find bin 
//                            node.GroupGain)

//    /// input data: groups of items
//    let dataGroups = dataGroupsA |> Array.map snd

//    /// initial configuration (after regular KMean)
//    let initialConfig = 
//        initConfig
//        |> Array.map 
//            (Array.map (fun (label,_) -> 
//                dataGroupsA 
//                |> Array.findIndex (fun (l,_) -> l=label) )
//            )

//    /// mapping between group pairs (and matrices) and pq
//    let pairArrayG = pairArrayFn dataGroups.Length

//    //let fComp = ClusterCheck.checkSgain gainFn (matrixSingletons)

//    /// adjustency matrix (initialize with initial configuration) 
//    let matrixA =
//        let m = JaggedArray.zeroCreate dataGroups.Length dataGroups.Length
//        pairArrayG
//        |> Array.iter (fun (i,j) -> 
//            let cluster = initialConfig |> Array.find (fun x -> x |> Array.contains i)
//            m.[i].[i] <- 1
//            if (cluster |> Array.contains j) then   
//                m.[i].[j] <- 1
//            else
//                m.[i].[j] <- 0
//            )
//        m
//        |> Array.map (Array.toList)
//        |> Array.toList

//    let rename (list: (string*(Types.Item [])) []) =
//            if list.Length>1 then
//                let newKey =
//                    list 
//                    |> Array.fold (fun state (k,v) -> (sprintf "%s|%s" state k)) (sprintf "mix")  
//                let newListValue = 
//                    list 
//                    |> Array.map (fun (k,v) -> v) 
//                    |> Array.concat      
//                    |> Array.map (fun protein -> {protein with BinL= Array.append protein.BinL.[0 .. depth] [|newKey|]})   
//                (newKey, newListValue)
//            else
//                list.[0]    

//    superFunctionTestG dataGroups singleGG gainFn matrixSingletons matrixA
//    |> Seq.maxBy (fun (x,_) -> x)
//    |> snd
//    |> Array.map (fun groupIDs ->     
//        groupIDs 
//        |> Array.map (fun groupID -> 
//            dataGroupsA.[groupID])
//        |> rename
//            )
//    |> Map.ofArray

////!!! Rethink idea of directions!!! not element -> element but element -> cluster

//let applySSNcombi data setN = createTreeWalking (SSN.getStepGainNodeSetnR setN) None SSN_combi data
//let applySSN data setN = createTreeWalking (SSN.getStepGainNodeSetnR setN) None (SSN_walk (FunctionsExp.PreClusterFunctions.kmeanGroupsKKZ, (walkingFn))) data
//let readMM data setN = createTreeWalking (SSN.getStepGainNodeSetnR setN) None MM data

//applySSN (TestData.SynData.data1) 10
//applySSNcombi (TestData.SynData.data1) 10

//stepCount

//let data6 = 
//    ArabiProteome.itemsWithMapManFound
//    |> Array.filter (fun x -> x.BinL.[0]="6")
//    |> Array.mapi (fun id x -> {x with ID=id})

//let treeMM = readMM data6 data6.Length
//treeMM.GroupGain
//let walk = applySSN (data6) data6.Length
//walk.GroupGain
//let combi = applySSNcombi (data6) data6.Length
//combi.GroupGain

//let matrixOfN = 
//    data6
//    |> distMatrixWeightedOf distanceMatrixWeighted None
//    /// normalised
//let dissMax = 
//    matrixOfN
//    |> Array2D.array2D_to_seq
//    |> Seq.max

//let normMatrix = 
//    matrixOfN
//    |> Array2D.map (fun i -> i/dissMax)

//walk
//|> Tree.filterLeaves
//|> Analysis.pointDxC normMatrix
//|> fun (x,y) -> weightedEuclidean None [x;y] [0.;0.]

//let timeAra = [|15.;180.;2880.;5760.;(5760.+15.);(5760.+180.);(5760.+2880.);(5760.+5760.)|]
//Plots.drawKinetik (walk |> Tree.findNode ["6"]) timeAra "bin 6" |> Chart.Show//AsImage (StyleParam.ImageFormat.JPEG)

//treeMM

//GePhi.sendToGephiFromTreeParam treeMM

//GePhi.sendToGephiFromTreeParam walk
//GePhi.sendToGephiFromTreeParam combi


//let pathsSorted =
//    ArabiProteome.itemsWithMapManFound
//    |> Array.groupBy (fun x -> x.BinL.[0])
//    |> Array.filter (fun (bin,l) -> l.Length > 2)
//    |> Array.sortBy (fun (bin,l) -> l.Length )
//    |> Array.map (fst)
    
//let pathsSortedChildren =
//    ChlamyProteome.dataAll //ArabiProteome.itemsWithMapManFound
//    |> Array.groupBy (fun x -> x.BinL.[0])
//    |> Array.filter (fun (bin,l) -> l.Length > 2)
//    |> Array.map (fun (bin,l) -> 
//        let data = l |> Array.mapi (fun id x -> {x with ID=id})
//        let tree = readMM data data.Length
//        let childrenN = Tree.filterChildrenList tree
//        (bin, childrenN |> List.max) )
//    |> Array.filter (fun (bin,cn) -> cn<=10)
//    |> Array.map fst

//let pathFile = @"c:\Users\mikha\Work-CSB\Projects\SSN\results\walkResults\"
//let neader = "Path\tcombi_GG\tcombi_Time\tcombi_DxC\twalk_GG\twalk_Time\twalk_DxC\tcomparison"

//let lines = 
//    [|"1";"31"|]
//    |> List.ofArray
//    |> List.map (fun path -> 
//        let data = 
//            ChlamyProteome.dataAll //ArabiProteome.itemsWithMapManFound
//            |> Array.filter (fun x -> x.BinL.[0]=path)
//            |> Array.mapi (fun id x -> {x with ID=id})
//        let stopwatch = new System.Diagnostics.Stopwatch()
//        stopwatch.Start()
//        let walk = applySSN (data) data.Length
//        let timeWalk = (stopwatch.Elapsed.TotalSeconds)
//        let combi = applySSNcombi (data) data.Length
//        let timeCombi = (stopwatch.Elapsed.TotalSeconds)
//        stopwatch.Stop()
//        let matrixOfN = 
//            data
//            |> distMatrixWeightedOf distanceMatrixWeighted None
//            /// normalised
//        let dissMax = 
//            matrixOfN
//            |> Array2D.array2D_to_seq
//            |> Seq.max

//        let normMatrix = 
//            matrixOfN
//            |> Array2D.map (fun i -> i/dissMax)

//        let dxcW = 
//            walk
//            |> Tree.filterLeaves
//            |> Analysis.pointDxC normMatrix
//            |> fun (x,y) -> weightedEuclidean None [x;y] [0.;0.]
    
//        let dxcC = 
//            combi
//            |> Tree.filterLeaves
//            |> Analysis.pointDxC normMatrix
//            |> fun (x,y) -> weightedEuclidean None [x;y] [0.;0.]

//        printfn "path: %s, %i" path data.Length 
//        printfn "combi GG; combi Time; walk GG; walk Time; walk DxC; comparison"
//        printfn "%s %f %f %f %f %f %f %i" path combi.GroupGain timeCombi dxcC walk.GroupGain timeWalk dxcW (Tree.treeComparison combi walk)
//        sprintf "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%i" path combi.GroupGain timeCombi dxcC walk.GroupGain timeWalk dxcW (Tree.treeComparison combi walk)

//        //(path, data.Length, combi.GroupGain, timeCombi, dxcC, walk.GroupGain, timeWalk, dxcW, (Tree.treeComparison combi walk))
//        )


//let linesWAcombi = 
//    [|"11"|]
//    |> List.ofArray
//    |> List.map (fun path -> 
//        let data = 
//            ArabiProteome.itemsWithMapManFound
//            |> Array.filter (fun x -> x.BinL.[0]=path)
//            |> Array.mapi (fun id x -> {x with ID=id})
//        let stopwatch = new System.Diagnostics.Stopwatch()
//        stopwatch.Start()
//        let walk = applySSN (data) data.Length
//        let timeWalk = (stopwatch.Elapsed.TotalSeconds)
//        //let combi = applySSNcombi (data) data.Length
//        //let timeCombi = (stopwatch.Elapsed.TotalSeconds)
//        stopwatch.Stop()
//        let matrixOfN = 
//            data
//            |> distMatrixWeightedOf distanceMatrixWeighted None
//            /// normalised
//        let dissMax = 
//            matrixOfN
//            |> Array2D.array2D_to_seq
//            |> Seq.max

//        let normMatrix = 
//            matrixOfN
//            |> Array2D.map (fun i -> i/dissMax)

//        let dxcW = 
//            walk
//            |> Tree.filterLeaves
//            |> Analysis.pointDxC normMatrix
//            |> fun (x,y) -> weightedEuclidean None [x;y] [0.;0.]
    
//        //let dxcC = 
//        //    combi
//        //    |> Tree.filterLeaves
//        //    |> Analysis.pointDxC normMatrix
//        //    |> fun (x,y) -> weightedEuclidean None [x;y] [0.;0.]

//        printfn "path: %s, %i" path data.Length 
//        printfn "combi GG; combi Time; walk GG; walk Time; walk DxC; comparison"
//        printfn "%s %f %f %f %f %f %f %i" path 0. 0. 0. walk.GroupGain timeWalk dxcW 0
//        sprintf "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%i" path 0. 0. 0. walk.GroupGain timeWalk dxcW 0

//        //(path, data.Length, combi.GroupGain, timeCombi, dxcC, walk.GroupGain, timeWalk, dxcW, (Tree.treeComparison combi walk))
//        )

//File.AppendAllLines((sprintf "%s%s.txt" pathFile "oldChlamyProtData"), neader :: lines)

