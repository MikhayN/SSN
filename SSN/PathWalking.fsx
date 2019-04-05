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

type ClusterFn = int -> Map<string, Item []> -> (string * (Item [])) [] []

type WalkingFn = (float -> float -> int -> int -> float) -> float [,] -> Map<string, Item []> -> (string * (Item [])) [] [] -> Map<string, Item []>

type Mode =
    |MM                                     // MapMan with broken leaves
    |SSN_combi                              // without simplification, pure combinatorics
    |SSN_walk of (ClusterFn*WalkingFn)      // with kMean as a start point for gain walking


let createTreeWalking gainFn (weight: seq<float> option) (mode: Mode) (rootGroup: Types.Item array) = 
        
    let nRoot = rootGroup.Length

    let matrixSingletons = 
        rootGroup
        |> distMatrixWeightedOf distanceMatrixWeighted weight

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
                        [|2 .. (nMax-1)|] 
                        |> Array.map (fun i -> 
                            x
                            |> kmeanKKZ i
                            |> walkFn gainFn matrixSingletons x) )
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

let pairwiseCorrAverage (x:matrix) (y:matrix) =
    let xN = x.Dimensions |> fst
    let yN = y.Dimensions |> fst
    let m = Array2D.create xN yN 0.
    for rowI in 0..(xN-1) do
        for colI in 0..(yN-1) do
            let tmp = Functions.SSN.weightedEuclidean None (x.Row rowI) (y.Row colI) 
            m.[rowI,colI] <- tmp
    m 
    |> Array2D.array2D_to_seq 
    |> Seq.average

let distMatrix (matrixA: string * matrix) (matrixB: string * matrix) =
    let mA = snd matrixA
    let mB = snd matrixB
    pairwiseCorrAverage mA mB

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

let updateCentroid (current: string * matrix) (sample: (string * matrix) []) = 
    let size = sample.Length
    match size with
    | 0 -> current
    | _ ->
        ("", 
            sample
            |> Array.map (fun (_,x) -> x.ToArray2D() |> Array2D.toJaggedArray)
            |> Array.concat
            |> fun p -> MatrixTopLevelOperators.matrix p)
    
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


/// New functions for PathWalk

let pairArrayFn (singletons: int []) = 
    singletons 
    |> Array.allPairs singletons
    |> Array.map (fun (a,b) -> if a>b then (b,a) else (a,b))
    |> Array.distinct

/// input: i - index of a foreign element, clusterBinaryVector - binary vector, representing a cluster.
let dist (i: int) (clusterBinaryVector: int []) (matrixD: float [,]) =
    //in clusterBinaryVector if 1 - elements are in the same cluster, 0 - foreigners)
    clusterBinaryVector 
    |> Array.indexed 
    |> Array.filter (fun (_, clusterBin) ->  clusterBin=1) 
    |> Array.averageBy (fun (j, _) -> matrixD.[i,j])


let reflectState (matrixA) (data: Item []) =
    matrixA 
    |> List.distinct
    |> List.map (fun i -> 
        i 
        |> List.indexed 
        |> List.filter (fun (_,x) -> x=1 ) 
        |> List.map (fun (i,_) -> data.[i])
        |> List.toArray
        )
    |> List.toArray
    
/// loop for path walking with nDirections iterations and nSteps in depth, applying f on each step
let superFunctionTest nDirections nSteps data matrixD (pairArray: (int*int) []) fComp matrixA pq_origin =

    let rec loop iStep (mA: int list list) (pq: MinIndexPriorityQueue<float>)  =
        
        seq [for iDirection in [1 .. nDirections] do
                
                printfn "direction %i at depth %i" iDirection iStep
                printfn "scheme before: %A" (data |> reflectState mA |> Array.map groupIDFn)

                let pairID = pq.TopIndex()
                pq.Pop() |> ignore          // pq will be used for other directiones in for loop
                let pq_new = pq.DeepCopy()  // pq_new will go for next step deeper in recursive part

                let (a,b) = // order represents the moving: a - moved, b - target
                        pairArray.[pairID] 
                        |> (fun (a,b) -> 
                            let B = mA.[b] |> List.toArray
                            let A = mA.[a] |> List.toArray
                            if (dist a B matrixD) < (dist b A matrixD) then  
                                a,b
                            else
                                b,a
                            )
                printfn "moving: %i -> %i" a b

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

                /// updatePQ
                let pairs_aB =  
                    B |> Array.map (fun i -> 
                        if a>i then (i,a) else (a,i)
                        |> fun (a,iB) -> pairArray |> Array.findIndex (fun pair -> pair=(a,iB))
                        )
                let pairs_aA = 
                    A 
                    |> Array.filter (fun i -> i<>a)
                    |> Array.map (fun i -> 
                        if a>i then (i,a) else (a,i)
                        |> fun (a,iA) -> pairArray |> Array.findIndex (fun pair -> pair=(a,iA))
                        )
                pq_new.TryRemoveGroup pairs_aB
                pq_new.TryReturnGroup pairs_aA

                let mA_new = matrixAA |> Array.map (Array.toList) |> Array.toList
                                
                let newState = reflectState mA_new data
                let gain = fComp (newState)
                let configuration = (newState) |> Array.map groupIDFn
                let stats = (gain, configuration, iDirection, iStep)
                
                printfn "%f, %A" gain configuration 

                if iStep=nSteps then
                    yield (stats)
                else
                    yield (stats)
                    yield! loop (iStep+1) mA_new pq_new

        ]
    
    let r = loop 1 matrixA pq_origin
        
    r |> Seq.toList

/// Test Path Walk

let singletons = TestData.SynData.data1

let data = singletons
//Plots.drawKinetik data [|0 .. 5|] "" |> Chart.Show
let kMeanScheme = [| [|data.[0]; data.[2]|] ; [|data.[1]; data.[3]|] ; [|data.[4]|] ; [|data.[5]; data.[7]; data.[8]|] ; [|data.[6]|] ; [|data.[9]|] |]

let matrixD = data |> distMatrixWeightedOf distanceMatrixWeighted None
let fComp = ClusterCheck.checkSgain (SSN.getStepGainNodeSetnR 10) (matrixD)

let pairArray = pairArrayFn (groupIDFn data)

let pq_o =
    let n = pairArray.Length
    let p = MinIndexPriorityQueue<float>(n) // n-i
    for j=0 to n-1 do 
        let (i,ii) = pairArray.[j]
        p.Insert j matrixD.[i,ii]
    p

let matrixA =
    let m = JaggedArray.zeroCreate data.Length data.Length
    pairArray
    |> Array.iteri (fun id (i,j) -> 
        let cluster = kMeanScheme |> Array.find (fun x -> x |> Array.contains (data.[i]))
        if (cluster |> Array.contains (data.[j])) then   
            printfn "remove id=%i pair %i %i" id i j
            pq_o.Remove id
            m.[i].[j] <- 1
            m.[j].[i] <- 1
        else
            m.[i].[j] <- 0
            m.[j].[i] <- 0 
        )
    m
    |> Array.map (Array.toList)
    |> Array.toList
    
let nDirections = 3
let nSteps = 6
let nCalc = [1 .. nSteps] |> List.sumBy (fun i -> (float nDirections)**(float i)) |> int    
let dataToPlot = (fComp kMeanScheme, kMeanScheme |> Array.map groupIDFn, 0 , 0) :: (superFunctionTest 3 5 data matrixD pairArray fComp matrixA pq_o)

dataToPlot.Length 
dataToPlot |> List.maxBy (fun (x,_,_,_) -> x)
dataToPlot |> List.filter (fun (x,_,_,_) -> x>=41.41360694)

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
                                    i - (list 
                                    |> List.cutAfterN (i+1)
                                    |> fst
                                    |> List.rev
                                    |> List.findIndex (fun (_,_,_,dBefore) -> dBefore < dCurrent)) 
                                let target = i+1
                                
                                (source,target))
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


plotGainWalk dataToPlot

//
//// Adjust the PathWalking to the groups instead of singletons
//

let distanceMatrixWeighted f (dataKinetics: float [] [] []) =
        let m = Array2D.zeroCreate (dataKinetics.Length) (dataKinetics.Length)
        for rowI in 0..dataKinetics.Length-1 do
            for colI in 0..rowI do
                let tmp = f (matrix dataKinetics.[rowI]) (matrix dataKinetics.[colI]) 
                m.[colI,rowI] <- tmp
                m.[rowI,colI] <- tmp
        m

let pairArrayFnG (groupIDs: int []) = 
    groupIDs 
    |> Array.allPairs groupIDs
    |> Array.map (fun (a,b) -> if a>b then (b,a) else (a,b))
    |> Array.distinct

/// input: i - index of a foreign element, clusterBinaryVector - binary vector, representing a cluster (1 inside, 0 outside).
let distG (i: int) (clusterBinaryVector: int []) (matrixG: float [,]) =
    //in clusterBinaryVector if 1 - elements are in the same cluster, 0 - foreigners)
    clusterBinaryVector 
    |> Array.indexed 
    |> Array.filter (fun (_, clusterBin) ->  clusterBin=1) 
    |> Array.averageBy (fun (j, _) -> matrixG.[i,j])


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
let superFunctionTestG nDirections nSteps data matrixD (pairArray: (int*int) []) fComp matrixA pq_origin =

    let rec loop iStep (mA: int list list) (pq: MinIndexPriorityQueue<float>)  =
        
        seq [for iDirection in [1 .. nDirections] do
                
                printfn "direction %i at depth %i" iDirection iStep
                printfn "scheme before: %A" (data |> reflectStateG mA |> Array.map groupIDFn)

                let pairID = pq.TopIndex()
                pq.Pop() |> ignore          // pq will be used for other directiones in for loop
                let pq_new = pq.DeepCopy()  // pq_new will go for next step deeper in recursive part

                let (a,b) = // order represents the moving: a - moved, b - target
                        pairArray.[pairID] 
                        |> (fun (a,b) -> 
                            let B = mA.[b] |> List.toArray
                            let A = mA.[a] |> List.toArray
                            if (distG a B matrixD) < (distG b A matrixD) then  
                                a,b
                            else
                                b,a
                            )
                printfn "moving: %i -> %i" a b

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

                /// updatePQ
                let pairs_aB =  
                    B |> Array.map (fun i -> 
                        if a>i then (i,a) else (a,i)
                        |> fun (a,iB) -> pairArray |> Array.findIndex (fun pair -> pair=(a,iB))
                        )
                let pairs_aA = 
                    A 
                    |> Array.filter (fun i -> i<>a)
                    |> Array.map (fun i -> 
                        if a>i then (i,a) else (a,i)
                        |> fun (a,iA) -> pairArray |> Array.findIndex (fun pair -> pair=(a,iA))
                        )
                pq_new.TryRemoveGroup pairs_aB
                pq_new.TryReturnGroup pairs_aA

                let mA_new = matrixAA |> Array.map (Array.toList) |> Array.toList
                                
                let newState = reflectStateG mA_new data
                let gain = fComp (newState)
                let configuration = (newState) |> Array.map groupIDFn
                let stats = (gain, configuration, iDirection, iStep)
                
                printfn "%f, %A" gain configuration 

                if iStep=nSteps then
                    yield (stats)
                else
                    yield (stats)
                    yield! loop (iStep+1) mA_new pq_new

        ]
    
    let r = loop 1 matrixA pq_origin
        
    r |> Seq.toList

let f = (SSN.getStepGainNodeSetnR 10)

/// solid groups of items
let dataGroups = singletons |> Array.map (fun i -> i.ProteinL.[0], [|i|])

//Plots.drawKinetik data [|0 .. 5|] "" |> Chart.Show
let kMeanSchemeInit = [| [|dataGroups.[0]; dataGroups.[2]|] ; [|dataGroups.[1]; dataGroups.[3]|] ; [|dataGroups.[4]|] ; [|dataGroups.[5]; dataGroups.[7]; dataGroups.[8]|] ; [|dataGroups.[6]|] ; [|dataGroups.[9]|] |]

//let parentGroup = kMeanSchemeInit |> Array.map (Array.map snd >> Array.concat) |> Array.concat

//let gainDiff listA listB : float =
//    let gFn (current: Types.Item []) = SSN.getStepGainFn f (SSN.groupIDFn current) (SSN.groupIDFn parentGroup) parentGroup.Length matrixD
//    (gFn (Array.append listA listB))

/// here G stands for group, but the value is still distance
let matrixG =
    distanceMatrixWeighted pairwiseCorrAverage (dataGroups |> Array.map (snd >> Array.map (fun i -> i.dataL))) 

//let fComp = ClusterCheck.checkSgain f (matrixD)

let pairArrayG = pairArrayFnG [|0 .. (dataGroups.Length-1)|]

let pq_oG =
    let n = pairArrayG.Length
    let p = MinIndexPriorityQueue<float>(n) // n-i
    for j=0 to n-1 do 
        let (i,ii) = pairArrayG.[j]
        p.Insert j matrixG.[i,ii]
    p

let matrixAG =
    let m = JaggedArray.zeroCreate data.Length data.Length
    pairArray
    |> Array.iteri (fun id (i,j) -> 
        let cluster = kMeanSchemeInit |> Array.find (fun x -> x |> Array.contains (dataGroups.[i]))
        if (cluster |> Array.contains (dataGroups.[j])) then   
            printfn "remove id=%i pair %i %i" id i j
            pq_oG.Remove id
            m.[i].[j] <- 1
            m.[j].[i] <- 1
        else
            m.[i].[j] <- 0
            m.[j].[i] <- 0 
        )
    m
    |> Array.map (Array.toList)
    |> Array.toList
    
//let nDirections = 3
//let nSteps = 5
//let nCalc = [1 .. nSteps] |> List.sumBy (fun i -> (float nDirections)**(float i)) |> int    
let dataToPlotG = 
    (fComp kMeanScheme, kMeanScheme |> Array.map groupIDFn, 0 , 0) 
        :: (superFunctionTestG 3 5 (dataGroups |> Array.map snd) matrixG pairArrayG fComp matrixAG pq_oG)


dataToPlotG.Length 
dataToPlotG |> List.maxBy (fun (x,_,_,_) -> x)
dataToPlotG |> List.filter (fun (x,_,_,_) -> x>=41.41360694)

plotGainWalk dataToPlotG

///// change matrixGroup to matrixGain (asymmetrical)

////// IC = -(nC/nR)*log2(nC/nR)
let f__ setNR dCurrSum dPredSum numberCurr numberRoot =
    let nC = float numberCurr
    let nR = float setNR
    let deltaDist = dPredSum - dCurrSum
    let deltaSpec = -(nC/nR)*log2(nC/nR)//+((nR-nC)/nR)*log2((nR-nC)/nR))
    if numberCurr=numberRoot then
        0.
    else
        deltaDist*deltaSpec

//
//// Adjust the PathWalking to using the Gain instead of Distance measure
//

let distanceMatrixWeighted f (dataKinetics: float [] [] []) =
        let m = Array2D.zeroCreate (dataKinetics.Length) (dataKinetics.Length)
        for rowI in 0..dataKinetics.Length-1 do
            for colI in 0..rowI do
                let tmp = f (matrix dataKinetics.[rowI]) (matrix dataKinetics.[colI]) 
                m.[colI,rowI] <- tmp
                m.[rowI,colI] <- tmp
        m

let pairArrayFnGain (groupIDs: int []) = 
    groupIDs 
    |> Array.allPairs groupIDs

/// input: i - index of a foreign element, clusterBinaryVector - binary vector, representing a cluster (1 inside, 0 outside).
let distG (i: int) (clusterBinaryVector: int []) (matrixG: float [,]) =
    //in clusterBinaryVector if 1 - elements are in the same cluster, 0 - foreigners)
    clusterBinaryVector 
    |> Array.indexed 
    |> Array.filter (fun (_, clusterBin) ->  clusterBin=1) 
    |> Array.averageBy (fun (j, _) -> matrixG.[i,j])


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
let superFunctionTestG nDirections nSteps data matrixD (pairArray: (int*int) []) fComp matrixA pq_origin =

    let rec loop iStep (mA: int list list) (pq: MinIndexPriorityQueue<float>)  =
        
        seq [for iDirection in [1 .. nDirections] do
                
                printfn "direction %i at depth %i" iDirection iStep
                printfn "scheme before: %A" (data |> reflectStateG mA |> Array.map groupIDFn)

                let pairID = pq.TopIndex()
                pq.Pop() |> ignore          // pq will be used for other directiones in for loop
                let pq_new = pq.DeepCopy()  // pq_new will go for next step deeper in recursive part

                let (a,b) = // order represents the moving: a - moved, b - target
                        pairArray.[pairID] 
                        |> (fun (a,b) -> 
                            let B = mA.[b] |> List.toArray
                            let A = mA.[a] |> List.toArray
                            if (distG a B matrixD) < (distG b A matrixD) then  
                                a,b
                            else
                                b,a
                            )
                printfn "moving: %i -> %i" a b

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

                /// updatePQ
                let pairs_aB =  
                    B |> Array.map (fun i -> 
                        if a>i then (i,a) else (a,i)
                        |> fun (a,iB) -> pairArray |> Array.findIndex (fun pair -> pair=(a,iB))
                        )
                let pairs_aA = 
                    A 
                    |> Array.filter (fun i -> i<>a)
                    |> Array.map (fun i -> 
                        if a>i then (i,a) else (a,i)
                        |> fun (a,iA) -> pairArray |> Array.findIndex (fun pair -> pair=(a,iA))
                        )
                pq_new.TryRemoveGroup pairs_aB
                pq_new.TryReturnGroup pairs_aA

                let mA_new = matrixAA |> Array.map (Array.toList) |> Array.toList
                                
                let newState = reflectStateG mA_new data
                let gain = fComp (newState)
                let configuration = (newState) |> Array.map groupIDFn
                let stats = (gain, configuration, iDirection, iStep)
                
                printfn "%f, %A" gain configuration 

                if iStep=nSteps then
                    yield (stats)
                else
                    yield (stats)
                    yield! loop (iStep+1) mA_new pq_new

        ]
    
    let r = loop 1 matrixA pq_origin
        
    r |> Seq.toList

let parentGroup = dataGroups |> Array.map (snd) |> Array.concat

let gainDiff listA listB : float =
    let gFn (current: Types.Item []) = SSN.getStepGainFn f (SSN.groupIDFn current) (SSN.groupIDFn parentGroup) parentGroup.Length matrixD
    (gFn (Array.append listA listB))

/// gains between groups of elements, if i moved to cluster with j
// matrixGain can be fixed throughout the whole algorithm or be updated based on the previous cluster configuration
let matrixGain = 
    distanceMatrixWeighted pairwiseCorrAverage (dataGroups |> Array.map (snd >> Array.map (fun i -> i.dataL))) 

//let fComp = ClusterCheck.checkSgain f (matrixD)

let pairArrayGain = pairArrayFnGain [|0 .. (dataGroups.Length-1)|]

let pq_oGain =
    let n = pairArrayGain.Length
    let p = MinIndexPriorityQueue<float>(n) // n-i
    for j=0 to n-1 do 
        let (i,ii) = pairArrayG.[j]
        p.Insert j matrixGain.[i,ii]
    p

let matrixAGain =
    let m = JaggedArray.zeroCreate data.Length data.Length
    pairArrayGain
    |> Array.iteri (fun id (i,j) -> 
        let cluster = kMeanSchemeInit |> Array.find (fun x -> x |> Array.contains (dataGroups.[i]))
        if (cluster |> Array.contains (dataGroups.[j])) then   
            printfn "remove id=%i pair %i %i" id i j
            pq_oGain.Remove id
            m.[i].[j] <- 1
            m.[j].[i] <- 1
        else
            m.[i].[j] <- 0
            m.[j].[i] <- 0 
        )
    m
    |> Array.map (Array.toList)
    |> Array.toList
    
let nDirections = 3
let nSteps = 10
let nCalc = [1 .. nSteps] |> List.sumBy (fun i -> (float nDirections)**(float i)) |> int    

let dataToPlotGain = 
    (fComp kMeanScheme, kMeanScheme |> Array.map groupIDFn, 0 , 0) 
        :: (superFunctionTestG 3 5 (dataGroups |> Array.map snd) matrixGain pairArrayGain fComp matrixAGain pq_oGain)


dataToPlotGain.Length 
dataToPlotGain |> List.maxBy (fun (x,_,_,_) -> x)
dataToPlotGain |> List.filter (fun (x,_,_,_) -> x>=41.41360694)

plotGainWalk dataToPlotGain


/// Test 

//let walkingFunction nDir nSt (fGain) matrixD (initialConf: (string * (Item [])) [] []) : Map<string, Item []> =    
    
//    let data = singletons
//    //Plots.drawKinetik data [|0 .. 5|] "" |> Chart.Show
//    let kMeanScheme = [| [|data.[0]; data.[2]|] ; [|data.[1]; data.[3]|] ; [|data.[4]|] ; [|data.[5]; data.[7]; data.[8]|] ; [|data.[6]|] ; [|data.[9]|] |]

//    let matrixD = data |> distMatrixWeightedOf distanceMatrixWeighted None
//    let fComp = ClusterCheck.checkSgain (SSN.getStepGainNodeSetnR 10) (matrixD)

//    let pairArray = pairArrayFn (groupIDFn data)

//    let pq_o =
//        let n = pairArray.Length
//        let p = MinIndexPriorityQueue<float>(n) // n-i
//        for j=0 to n-1 do 
//            let (i,ii) = pairArray.[j]
//            p.Insert j matrixD.[i,ii]
//        p

//    let matrixA =
//        let m = JaggedArray.zeroCreate data.Length data.Length
//        pairArray
//        |> Array.iteri (fun id (i,j) -> 
//            let cluster = kMeanScheme |> Array.find (fun x -> x |> Array.contains (data.[i]))
//            if (cluster |> Array.contains (data.[j])) then   
//                printfn "remove id=%i pair %i %i" id i j
//                pq_o.Remove id
//                m.[i].[j] <- 1
//                m.[j].[i] <- 1
//            else
//                m.[i].[j] <- 0
//                m.[j].[i] <- 0 
//            )
//        m
//        |> Array.map (Array.toList)
//        |> Array.toList

    
//    let fComp = ClusterCheck.checkSgain fGain (matrixD)






//////// Usage
//let applySSNnew data setN = createTreeWalking (SSN.getStepGainNodeSetnR setN) None (SSN (kmeanGroupsKKZ, walkingFunction nDirections nSteps)) data