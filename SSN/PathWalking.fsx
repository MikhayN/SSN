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

type ClusterFn = int -> seq<float> option -> Map<string, Item []> -> Item [] [] []

open BioFSharp.Elements

type WalkingFn = (float -> float -> int -> int -> float) -> float [,] ->  Item [] [] [] -> Map<string, Item []>

type Mode =
    |MM                                     // MapMan with broken leaves
    |SSN of (ClusterFn*WalkingFn)                     // with KMean Swap (scheme-wise approximation)
    |SSN_combi                              // without simplification, pure combinatorics


let createTree gainFn (weight: seq<float> option) (mode: Mode) (rootGroup: Types.Item array) = 
        
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

            |SSN (kmeanKKZ, walkFn) -> // with simplification overall
                if (nodeMembers.Length=1) 
                        || (nodeMembers.Length=0) 
                        || (String.contains "|" (String.Concat nodeMembers.[0].BinL)) then 
                    Map.empty
                else 
                    (breakGroup nodeMembers depth) // get a map with original grouping as (subbin, Item)
                    |> (fun x -> 
                        let nMax = x.Count
                        [|0 .. (nMax-2)|] 
                        |> Array.map (fun i -> 
                            x
                            |> kmeanKKZ (nMax-i) weight
                            |> walkFn gainFn matrix)
                    )
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
    
let kmeanGroupsKKZ (k: int) weight (children: Map<string,Types.Item []> ) =

    let data = children |> Map.toArray |> Array.map (fun (s,ar) -> (s,ar |> Array.map (fun p -> p.dataL) |> MatrixTopLevelOperators.matrix))

    let clusters = ML.Unsupervised.IterativeClustering.compute distMatrix (initCgroupsKKZ) updateCentroid data k

    data
    |> Array.map (fun list -> (clusters.Classifier list |> fst),list )
    |> Array.groupBy (fst)
    |> Array.map (fun (_,list) -> 
        list 
        |> Array.map (fun (_,(bin,_)) -> (bin, children |> Map.find bin))
        )

let singletons = TestData.SynData.data1

let pairArrayFn (singletons: int []) = 
    singletons 
    |> Array.allPairs singletons
    |> Array.map (fun (a,b) -> if a>b then (b,a) else (a,b))
    |> Array.distinct

  

let dist (i: int) (clusterBinaryVector: int []) (matrixD: float [,]) =
    //in clusterBinaryVector if 1 - elements are in the same cluster, 0 - foreigners)
    clusterBinaryVector |> Array.averageBy (fun j -> matrixD.[i,j])

let updateMatrices (ab: int*int) (pairArray: (int*int) []) (matrixA: int [] []) (pq: MinIndexPriorityQueue<float>) =
    
    /// update MatrixA, so that a is moved from A cluster to b in B cluster 

    let (a,b) = ab
        
    let (A,B) = (
        (matrixA.[a] |> Array.indexed |> Array.filter (fun (i,x) -> x=1) |> Array.map fst), // a in A
        (matrixA.[b] |> Array.indexed |> Array.filter (fun (i,x) -> x=1) |> Array.map fst)) // b in B
    
    // move 'a' out of 'A', but leave it stay with itself
    for i in A do   
        matrixA.[a].[i] <- 0
        matrixA.[i].[a] <- 0
    matrixA.[a].[a] <- 1
    
    // move a in B
    for j in B do
        matrixA.[a].[j] <- 1
        matrixA.[j].[a] <- 1
    
    /// updatePQ
    let pairs_aB =  
        B |> Array.map (fun i -> 
            if a>i then (i,a) else (a,i)
            |> fun (a,iB) -> pairArray |> Array.findIndex (fun pair -> pair=(a,iB)))
    let pairs_aA = 
        A 
        |> Array.filter (fun i -> i<>a)
        |> Array.map (fun i -> 
            if a>i then (i,a) else (a,i)
            |> fun (a,iA) -> pairArray |> Array.findIndex (fun pair -> pair=(a,iA)))
    pq.TryRemoveGroup pairs_aB
    pq.TryReturnGroup pairs_aA

let reflectState (matrixA) (data: Item []) =
    matrixA 
    |> Array.distinct
    |> Array.map (fun i -> i |> Array.indexed |> Array.filter (fun (i, x) -> x=1 ) |> Array.map (fun (i,_) -> data.[i]))

let mainPipeline data (matrixA: int [][]) (matrixD: float [,]) (pairArray: (int*int) []) (pq: MinIndexPriorityQueue<float>) =
    let pairID = pq.TopIndex()
    pq.Pop() |> ignore          // pq will be used for other directiones in for loop
    let pq_new = pq.DeepCopy()  // pq_new will go for next step deeper in recursive part

    //printfn "before update pq: %i" pq.Length
    //printfn "before update pq_new: %i" pq_new.Length 

    let (a,b) =
        pairArray.[pairID] 
        |> (fun (a,b) -> 
            let B = matrixA.[b]
            let A = matrixA.[a]
            if (dist a B matrixD) < (dist b A matrixD) then  
                a,b
            else
                b,a
            )
    
    printfn "item to move: %i to item; %i" a b

    updateMatrices (a,b) pairArray matrixA pq_new

    //printfn "after update pq: %i" pq.Length
    //printfn "after update pq_new: %i" pq_new.Length 

    let state = reflectState matrixA data
    (state, pq_new)


let data = singletons
let dataMap = data |> Array.map (fun p -> p.ProteinL.[0],[|p|]) |> Map.ofArray
let pairArray = pairArrayFn (groupIDFn data)
let matrixD = data |> distMatrixWeightedOf distanceMatrixWeighted None
let kMeanScheme = kmeanGroupsKKZ 6 None (dataMap) |> Array.map (Array.map (fun (_,i) -> i) >> Array.concat)

let pq_origin =
    let n = pairArray.Length
    let p = MinIndexPriorityQueue<float>(n) // n-i
    for j=0 to n-1 do 
        let (i,ii) = pairArray.[j]
        p.Insert j matrixD.[i,ii]
    p

matrixD |> Array2D.array2D_to_seq |> Seq.min

pq_origin.Top()

pairArray.Length

pq_origin.Pop()
pq_origin.Length
[1 .. pq_origin.Length] |> List.map (fun i -> pq_origin.HeapItem i)

let matrixA =
    let m = JaggedArray.zeroCreate data.Length data.Length
    pairArray
    |> Array.iteri (fun id (i,j) -> 
        let cluster = kMeanScheme |> Array.find (fun x -> x |> Array.contains (data.[i]))
        if (cluster |> Array.contains (data.[j])) then   
            pq_origin.Remove id
            m.[i].[j] <- 1
            m.[j].[i] <- 1
        else
            m.[i].[j] <- 0
            m.[j].[i] <- 0
        )
    m

let a = 8
let b = 9

(matrixA.[a] |> Array.indexed |> Array.filter (fun (i,x) -> x=1) |> Array.map fst), // a in A
(matrixA.[b] |> Array.indexed |> Array.filter (fun (i,x) -> x=1) |> Array.map fst)) // b in B


///     
let fStep : (MinIndexPriorityQueue<float> -> ((Item [] []) * MinIndexPriorityQueue<float>)) = mainPipeline data matrixA matrixD pairArray

let fComp = ClusterCheck.checkSgain (SSN.getStepGainNodeSetnR 10) (matrixD)

fComp kMeanScheme

/// loop for path walking with nDirections iterations and nSteps in depth, applying f on each step
let superFunction nDirections nSteps fStep fComp pq =

    let rec loop d statepq =
        seq [for a in [1 .. nDirections] do
                let newStatepq = fStep (snd statepq)
                if d=nSteps then
                    yield (newStatepq)
                else
                    yield (newStatepq)
                    yield! loop (d+1) newStatepq
            ]
    
    let r = loop 1 ([||],pq)

    let nCalc = [1 .. nSteps] |> List.sumBy (fun i -> (float nDirections)**(float i)) |> int
    printfn "number of checked states: %i" nCalc

    r |> Seq.fold (fun bestState currentState -> 
        if (fComp (fst currentState)) > (fComp (fst bestState)) then currentState else bestState ) (Seq.head r)


superFunction 3 3 fStep fComp pq_origin |> fst |> (fun x -> x, fComp x)

/// loop for path walking with nDirections iterations and nSteps in depth, applying f on each step
let superFunctionTest nDirections nSteps fStep fComp (pq :  MinIndexPriorityQueue<float>) =

    let rec loop d statepq =
        seq [for a in [1 .. nDirections] do
                let newStatepq : (Item [] []) * MinIndexPriorityQueue<float> = fStep (snd statepq)
                printfn "direction %i at depth %i" a d
                printfn "%f, %A, %i" (fComp (fst newStatepq)) (newStatepq |> fst |> Array.map groupIDFn)  ((snd newStatepq).Length)
                if d=nSteps then
                    yield (newStatepq)
                else
                    yield (newStatepq)
                    yield! loop (d+1) newStatepq
            ]
    
    let r = loop 1 ([||],pq)

    let nCalc = [1 .. nSteps] |> List.sumBy (fun i -> (float nDirections)**(float i)) |> int
    printfn "number of checked states: %i" nCalc
    
    r 
    |> Seq.map (fun (currentState, pq_c) -> 
        (fComp (currentState)), currentState |> Array.map groupIDFn , pq_c.Length)
    |> Seq.toArray

superFunctionTest 3 3 fStep fComp pq_origin
    