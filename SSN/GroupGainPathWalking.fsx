#r @"c:\Users\mikha\source\repos\mathnet-numerics\src\Numerics\bin\Debug\netstandard2.0\MathNet.Numerics.dll"
#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\netstandard.dll"

#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\BioFSharp.dll"
#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\BioFSharp.IO.dll"
#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharpAux.dll"
#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharpAux.IO.dll"
#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\Newtonsoft.Json.dll"
#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharp.Stats.dll"
#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharp.Plotly.dll"
#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharpGephiStreamer.dll"

#load "Types.fs"
#load "PQ.fs"
#load "Functions.fs"
#load "FunctionsExp.fs"
#load "GePhi.fs"
#load "TestData.fs"
#load "SOM.fs"
#load "Auxilliary.fs"
#load "Plots.fs"

open System 
open FSharpAux
open FSharp.Plotly
open FSharp.Stats

open Functions.SSN
open TestData
open GePhi
open Types
open PQ
open Auxilliary
open Plots
open FunctionsExp

#time

//let setN = 10

//let gainFn = (SSN.getStepGainNodeSetnR setN)

//// for matrixG recalculation 2)

//type GainComponents = {
//    ParentDiss : float
//    CurrentDiss : float
//    IC : float
//    }

//let ic nC nR = -(nC/nR)*log2(nC/nR)+((nR-nC)/nR)*log2((nR-nC)/nR)
//let parentDist i nR matrixD = SSN.dSumFn [|i|] [|0 .. (nR-1)|] matrixD
//let currentDist i clusterj matrixD = SSN.dSumFn [|i|] clusterj matrixD

//let gainValueFn (gainA: GainComponents []) =
//    gainA
//    |> Array.sumBy (fun x -> (x.ParentDiss - x.CurrentDiss) * x.IC )

//let gainComp i (matrixA: int list list)  matrixD =
//    let nR = matrixA.Length 
//    let clusterj = matrixA.[i] |> List.indexed |> List.filter (fun (i, x) -> x=1) |> List.map fst |> Array.ofList
//    {
//    ParentDiss = parentDist i nR matrixD;
//    CurrentDiss = currentDist i clusterj matrixD;
//    IC = ic (float clusterj.Length) (float nR)
//    }

//// utils

//let pairwiseCorrAverage (x:matrix) (y:matrix) =
//    let xN = x.Dimensions |> fst
//    let yN = y.Dimensions |> fst
//    let m = Array2D.create xN yN 0.
//    for rowI in 0..(xN-1) do
//        for colI in 0..(yN-1) do
//            let tmp = Functions.SSN.weightedEuclidean None (x.Row rowI) (y.Row colI) 
//            m.[rowI,colI] <- tmp
//    m 
//    |> Array2D.array2D_to_seq 
//    |> Seq.average


//let distanceMatrixWeighted f (dataKinetics: float [] [] []) =
//        let m = Array2D.zeroCreate (dataKinetics.Length) (dataKinetics.Length)
//        for rowI in 0..dataKinetics.Length-1 do
//            for colI in 0..rowI do
//                let tmp = f (matrix dataKinetics.[rowI]) (matrix dataKinetics.[colI]) 
//                m.[colI,rowI] <- tmp
//                m.[rowI,colI] <- tmp
//        m

let pairArrayFn (groupL: int) = 



    [|0 .. (groupL-1)|] 
    |> Array.allPairs [|0 .. (groupL-1)|] 
    |> Array.filter (fun (i,j) -> i<>j)

///// input: i - index of a foreign element, clusterBinaryVector - binary vector, representing a cluster.
//let dist (i: int) (clusterBinaryVector: int []) (matrixD: float [,]) =
//    //in clusterBinaryVector if 1 - elements are in the same cluster, 0 - foreigners)
//    clusterBinaryVector 
//    |> Array.indexed 
//    |> Array.filter (fun (_, clusterBin) ->  clusterBin=1) 
//    |> Array.averageBy (fun (j, _) -> matrixD.[i,j])

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
    

let reflectStateG (matrixA) (data: Item [] []) =
    matrixA 
    |> List.distinct
    |> List.map (fun i -> 
        i 
        |> List.indexed 
        |> List.filter (fun (_,x) -> x=1 ) 
        |> List.map (fun (i,_) -> data.[i])
        |> List.toArray
        |> Array.concat
        )
    |> List.toArray
    

/// loop for path walking with nDirections iterations and nSteps in depth, applying f on each step
let superFunctionTestG nDirections nSteps (data: Item [] []) (pairArray: (int*int) []) fComp matrixA  =

    let matrixG_origin =
        let m = JaggedArray.zeroCreate data.Length data.Length
        pairArray
        |> Array.iter (fun (i,j) ->     
            // 1 approach: direct calculation
            let matrixAA = matrixA |> List.map (List.toArray) |> List.toArray

            let (A,B) = (
                    (matrixAA.[i] |> Array.indexed |> Array.filter (fun (i,x) -> x=1) |> Array.map fst), // a was in A
                    (matrixAA.[j] |> Array.indexed |> Array.filter (fun (i,x) -> x=1) |> Array.map fst)) // b was in B

            for ii in A do   // move 'a' out of 'A', but leave it stay with itself
                    matrixAA.[i].[ii] <- 0
                    matrixAA.[ii].[i] <- 0
            matrixAA.[i].[i] <- 1

            for jj in B do   // move a in B
                    matrixAA.[i].[jj] <- 1
                    matrixAA.[jj].[i] <- 1
            let newMA = matrixAA |> Array.map (Array.toList) |> Array.toList
                        
            let temp = (reflectStateG newMA data)

            m.[i].[j] <- fComp temp
            // 2 approach: via GainComponents
            )
        m

    let pq_origin =
        let n = pairArray.Length
        let p = MaxIndexPriorityQueue<float>(n)
        for j=0 to n-1 do 
            let (i,ii) = pairArray.[j]
            p.Insert j matrixG_origin.[i].[ii] // load all gains except of (a -> a)
        pairArray |> Array.iteri (fun id (i,j) -> if  matrixA.[i].[j]=1 then p.Remove id) 
        p

    let rec loop iStep (mA: int list list) (pq: MaxIndexPriorityQueue<float>)  =
        
        seq [for iDirection in [1 .. nDirections] do
                
                //printfn "direction %i at depth %i" iDirection iStep
                //printfn "scheme before: %A" (data |> reflectStateG mA |> Array.map groupIDFn)

                let pairID = pq.TopIndex()
                pq.Pop() |> ignore          // pq will be used for other directiones in for loop
                //let pq_new = pq.DeepCopy()  // pq_new will go for next step deeper in recursive part

                let (a,b) = // order represents the moving: a - moved, b - target
                        pairArray.[pairID] 
                        
                //printfn "moving: %i -> %i" a b

                let matrixAA = mA |> List.map (List.toArray) |> List.toArray

                let (A,B) = (
                        (matrixAA.[a] |> Array.indexed |> Array.filter (fun (i,x) -> x=1) |> Array.map fst), // a was in A
                        (matrixAA.[b] |> Array.indexed |> Array.filter (fun (i,x) -> x=1) |> Array.map fst)) // b was in B

                for i in A do   // move 'a' out of 'A', but leave it stay with itself
                        matrixAA.[a].[i] <- 0
                        matrixAA.[i].[a] <- 0
                matrixAA.[a].[a] <- 1

                for j in B do   // move a in B
                        matrixAA.[a].[j] <- 1
                        matrixAA.[j].[a] <- 1

                ///// updatePQ
                //let pairs_aB =  
                //    B |> Array.map (fun i -> 
                //        if a>i then (i,a) else (a,i)
                //        |> fun (a,iB) -> pairArray |> Array.findIndex (fun pair -> pair=(a,iB))
                //        )
                //let pairs_aA = 
                //    A 
                //    |> Array.filter (fun i -> i<>a)
                //    |> Array.map (fun i -> 
                //        if a>i then (i,a) else (a,i)
                //        |> fun (a,iA) -> pairArray |> Array.findIndex (fun pair -> pair=(a,iA))
                //        )
                //pq_new.TryRemoveGroup pairs_aB
                //pq_new.TryReturnGroup pairs_aA

                let mA_new = matrixAA |> Array.map (Array.toList) |> Array.toList
                    
                let matrixG =
                    let m = JaggedArray.zeroCreate data.Length data.Length
                    pairArray
                    |> Array.iteri (fun id (i,j) ->     
                        // 1 approach: direct calculation
                        let matrixAA = mA_new |> List.map (List.toArray) |> List.toArray

                        let (A,B) = (
                                (matrixAA.[i] |> Array.indexed |> Array.filter (fun (i,x) -> x=1) |> Array.map fst), // a was in A
                                (matrixAA.[j] |> Array.indexed |> Array.filter (fun (i,x) -> x=1) |> Array.map fst)) // b was in B

                        for ii in A do   // move 'a' out of 'A', but leave it stay with itself
                                matrixAA.[i].[ii] <- 0
                                matrixAA.[ii].[i] <- 0
                        matrixAA.[i].[i] <- 1

                        for jj in B do   // move a in B
                                matrixAA.[i].[jj] <- 1
                                matrixAA.[jj].[i] <- 1
                        let newMA = matrixAA |> Array.map (Array.toList) |> Array.toList

                        m.[i].[j] <- fComp (reflectStateG newMA data)
                        // 2 approach: via GainComponents
                        )
                    m    

                let pq_new = 
                    let n = pairArray.Length
                    let p = MaxIndexPriorityQueue<float>(n)
                    for j=0 to n-1 do 
                        let (i,ii) = pairArray.[j]
                        p.Insert j matrixG.[i].[ii] // load all gains except of (a -> a)
                    pairArray |> Array.iteri (fun id (i,j) -> if  mA_new.[i].[j]=1 then p.Remove id) 
                    p


                let newState = reflectStateG mA_new data
                let gain = fComp (newState)
                let configuration = reflectStateID mA_new
                let stats = (gain, configuration) //, iDirection, iStep)
                
                //printfn "%f, %A" gain configuration 

                if iStep=nSteps then
                    yield (stats)
                else
                    yield (stats)
                    yield! loop (iStep+1) mA_new pq_new

        ]
    
    let r = loop 1 matrixA pq_origin
        
    r //|> Seq.toList

///// input data: groups of items
//let dataGroups = TestData.SynData.data1 |> Array.map (fun i ->  [|i|]) // i.ProteinL.[0],

///// initial configuration (after regular KMean)
//let initialConfig = 
//    [|
//    [|dataGroups.[0]; dataGroups.[2]|] ;
//    [|dataGroups.[0]; dataGroups.[2]|] ;
//    |]
//    //[| 
//    //[|dataGroups.[0]; dataGroups.[2]|] ; 
//    //[|dataGroups.[1]; dataGroups.[3]|] ; 
//    //[|dataGroups.[4]|] ; 
//    //[|dataGroups.[5]; dataGroups.[7]; dataGroups.[8]|] ; 
//    //[|dataGroups.[6]|] ; 
//    //[|dataGroups.[9]|] 
//    //|]

///// mapping between group pairs (and matrices) and pq
//let pairArrayG = pairArrayFn dataGroups.Length

///// intergroup distance matrix
//let matrixD =
//    distanceMatrixWeighted pairwiseCorrAverage (dataGroups |> Array.map (Array.map (fun i -> i.dataL))) 

//let fComp = ClusterCheck.checkSgain gainFn (matrixD)

///// adjustency matrix (initialize with initial configuration) 
//let matrixA =
//    let m = JaggedArray.zeroCreate dataGroups.Length dataGroups.Length
//    pairArrayG
//    |> Array.iteri (fun id (i,j) -> 
//        let cluster = initialConfig |> Array.find (fun x -> x |> Array.contains (dataGroups.[i]))
//        m.[i].[i] <- 1
//        if (cluster |> Array.contains (dataGroups.[j])) then   
//            m.[i].[j] <- 1
//        else
//            m.[i].[j] <- 0
//        )
//    m
//    |> Array.map (Array.toList)
//    |> Array.toList

///// matrix of gains after moving i -> j
//let matrixG =
//    let m = JaggedArray.zeroCreate dataGroups.Length dataGroups.Length
//    pairArrayG
//    |> Array.iteri (fun id (i,j) ->     
//        // 1 approach: direct calculation
//        let matrixAA = matrixA |> List.map (List.toArray) |> List.toArray

//        let (A,B) = (
//                (matrixAA.[i] |> Array.indexed |> Array.filter (fun (i,x) -> x=1) |> Array.map fst), // a was in A
//                (matrixAA.[j] |> Array.indexed |> Array.filter (fun (i,x) -> x=1) |> Array.map fst)) // b was in B

//        for ii in A do   // move 'a' out of 'A', but leave it stay with itself
//                matrixAA.[i].[ii] <- 0
//                matrixAA.[ii].[i] <- 0
//        matrixAA.[i].[i] <- 1

//        for jj in B do   // move a in B
//                matrixAA.[i].[jj] <- 1
//                matrixAA.[jj].[i] <- 1
//        let newMA = matrixAA |> Array.map (Array.toList) |> Array.toList

//        m.[i].[j] <- fComp (reflectStateG newMA dataGroups)
//        // 2 approach: via GainComponents
//        )
//    m

//let pq =
//    let n = pairArrayG.Length
//    let p = MaxIndexPriorityQueue<float>(n)
//    for j=0 to n-1 do 
//        let (i,ii) = pairArrayG.[j]
//        p.Insert j matrixG.[i].[ii] // load all gains except of (a -> a)
//    pairArrayG |> Array.iteri (fun id (i,j) -> if  matrixA.[i].[j]=1 then p.Remove id) 
//    p

//pq.Top()
//pairArrayG.[pq.TopIndex()]
//// try walk

//let f = (SSN.getStepGainNodeSetnR 10)

//let nDirections = 2
//let nSteps = 10
//let nCalc = [1 .. nSteps] |> List.sumBy (fun i -> (float nDirections)**(float i)) |> int    
//let kMeanScheme = initialConfig |> Array.map (Array.concat)

//let dataToPlot = 
//    (fComp kMeanScheme, kMeanScheme |> Array.map groupIDFn, 0 , 0) 
//    :: (superFunctionTestG nDirections nSteps dataGroups pairArrayG fComp matrixA)

//dataToPlot.Length 
//dataToPlot |> List.maxBy (fun (x,_,_,_) -> x)
//dataToPlot |> List.filter (fun (x,_,_,_) -> x>=41.41360694)


let plotGainWalk (list: ('a * 'b * int * int) list ) =
    let color = function
        |1 -> colorBlue
        |2 -> colorBrightBlue
        |3 -> colorGreen
        |4 -> colorOrange
        |5 -> colorYellow
        |_ -> colorGray
    
    let edges = 
        list 
        |> List.tail
        |> List.mapi ( fun i (_,_,_,dCurrent) -> 
                                
                                let source = 
                                    (list 
                                    |> List.cutAfterN (i+1)
                                    |> fst
                                    |> List.rev
                                    |> List.tryFindIndex (fun (_,_,_,dBefore) -> dBefore < dCurrent)) 
                                let target = i+1
                                match source with
                                |None -> (0,target)
                                |Some x -> (i - x,target)
                                )
    let lines = 
        edges
        |> List.map (fun (i,j) -> 
                        let (gain1,_,_,step1) = list.[i]
                        let (gain2,_,_,step2) = list.[j]
                        Chart.Line ([(step1,gain1);(step2,gain2)], Color = colorGray)
                    )
    let points =
        list
        |> List.map (fun (gain,_,direction,step) -> 
                        Chart.Point ([(step,gain)], Labels = [sprintf "%i" direction] ) |> Chart.withMarkerStyle (Size=10, Color=color direction)
                        )

    [lines; points] |> List.concat |> Chart.Combine |> Chart.Show

//plotGainWalk dataToPlot



// incorporate the pathWalking into the createTree


type ClusterFn = int -> Map<string, Item []> -> (string * (Item [])) [] []

type WalkingFn = int -> float [,] -> (float -> float -> int -> int -> float) -> Map<string, Item []> -> (string * (Item [])) [] [] -> Map<string, Item []>

type Mode =
    |MM                                     // MapMan with broken leaves
    |SSN_combi                              // without simplification, pure combinatorics
    |SSN_walk of (ClusterFn*WalkingFn)      // with kMean as a start point for gain walking


let createTreeWalking gainFn (weight: seq<float> option) (mode: Mode) (rootGroup: Types.Item array) = 
        
    let nRoot = rootGroup.Length

    let matrixSingletons = 
        rootGroup
        |> distMatrixWeightedOf SSN.distanceMatrixWeighted weight

    // calculation for one node    
    let rec loop (nodeMembers: Types.Item array) depth dPredSum =
        
        /// sum of max dist within a node 
        let dCurrSum = dSumFn (groupIDFn nodeMembers) (groupIDFn nodeMembers) matrixSingletons

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
                                    let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrixSingletons
                                    (loop nodes (depth+1) dPredSum'))


            |SSN_combi -> // without simplification, pure combinatorics
                if (nodeMembers.Length=1) 
                        || (nodeMembers.Length=0) 
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
                                    let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrixSingletons
                                    state 
                                    |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)
                            else
                                i
                                |> Map.fold (fun state key nodes ->
                                    match (singles.TryFind key) with
                                    | None ->
                                        let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrixSingletons
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

            |SSN_walk (kmeanKKZ, walkFn) -> // with kMean as a start point for gain walking
                if (nodeMembers.Length=1) 
                        || (nodeMembers.Length=0) 
                        || (String.contains "|" (String.Concat nodeMembers.[0].BinL)) then 
                    Map.empty
                else 
                    (breakGroup nodeMembers depth)
                    |> (fun x -> 
                        let nMax = x.Count
                        //printfn "length %i" nMax 
                        x::(
                            [2 .. (nMax-1)] 
                            |> List.rev
                            |> List.map (fun i -> 
                                x
                                |> kmeanKKZ i
                                |> walkFn depth matrixSingletons gainFn x) ) )// here the first in seq should be x, then everything else
                    |> Seq.fold (fun (singles,best) i -> 
                        let newNodes = 
                            if singles=Map.empty then
                                i
                                |> Map.fold (fun state key nodes ->
                                    let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrixSingletons
                                    state 
                                    |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)
                            else
                                i
                                |> Map.fold (fun state key nodes ->
                                    match (singles.TryFind key) with
                                    | None ->
                                        let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrixSingletons
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

let walkingFn nDirections nSteps depth matrixSingletons gainFn (data: Map<string, Item []>) (initConfig: (string * (Item [])) [] [])  = 
    
    let dataGroupsA = data |> Map.toArray

    /// input data: groups of items
    let dataGroups = dataGroupsA |> Array.map snd

    /// initial configuration (after regular KMean)
    let initialConfig = 
        initConfig
        |> Array.map 
            (Array.map (fun (label,_) -> 
                dataGroupsA 
                |> Array.findIndex (fun (l,_) -> l=label) )
            )

    /// mapping between group pairs (and matrices) and pq
    let pairArrayG = pairArrayFn dataGroups.Length

    ///// intergroup distance matrix
    //let matrixD =
    //    distanceMatrixWeighted pairwiseCorrAverage (dataGroups |> Array.map (Array.map (fun i -> i.dataL))) 

    //printfn "matrixD: %f" (sqrt (float matrixD.Length))

    let fComp = ClusterCheck.checkSgain gainFn (matrixSingletons)

    /// adjustency matrix (initialize with initial configuration) 
    let matrixA =
        let m = JaggedArray.zeroCreate dataGroups.Length dataGroups.Length
        pairArrayG
        |> Array.iteri (fun id (i,j) -> 
            let cluster = initialConfig |> Array.find (fun x -> x |> Array.contains i)
            m.[i].[i] <- 1
            if (cluster |> Array.contains j) then   
                m.[i].[j] <- 1
            else
                m.[i].[j] <- 0
            )
        m
        |> Array.map (Array.toList)
        |> Array.toList

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

    superFunctionTestG nDirections nSteps dataGroups pairArrayG fComp matrixA
    |> Seq.maxBy (fun (x,_) -> x)
    |> snd
    |> Array.map (fun groupIDs ->     
        groupIDs 
        |> Array.map (fun groupID -> 
            dataGroupsA.[groupID])
        |> rename
            )
    |> Map.ofArray

let applySSNcombi data setN = createTreeWalking (SSN.getStepGainNodeSetnR setN)  (None) SSN_combi data
let applySSN data setN = createTreeWalking (SSN.getStepGainNodeSetnR setN) (None) (SSN_walk (FunctionsExp.PreClusterFunctions.kmeanGroupsKKZ, (walkingFn 5 5))) data

applySSN (TestData.SynData.data4) 10
applySSNcombi (TestData.SynData.data4) 10

let data13 = 
    ChlamyProteome.dataAll
    |> Array.filter (fun x -> x.BinL.[0]="13")
    |> Array.mapi (fun id x -> {x with ID=id})

let walk = applySSN (data13) data13.Length
let combi = applySSNcombi (data13) data13.Length


Plots.drawKinetik (combi |> Tree.findNode ["1";"1";"mix|3|40"]) [|1.;24.;25.;26.;28.;32.|]  "subbin 1.1.mix|3|40" |> Chart.Show

GePhi.sendToGephiFromTreeParam walk
GePhi.sendToGephiFromTreeParam combi
