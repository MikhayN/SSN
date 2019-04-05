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
open System.IO 
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

/// file subname for the trial # change each time for the new file!!!!! ##############################################################################
//let fileSubName = "ProteinDataSet_Path8_noRep__nostates_plusInit_d3s7"

///// file for tracking a path during the whole treeCreate. 
//let header = "StepN\tDirectionN\tmovedID\tMovedTo\tGainBefore\tGainNow" 
//let pathLOG = @"c:\Users\mikha\Work-CSB\Projects\SSN\results\walkResults\PathLog_"

//File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [header])

///
let fileLogName = "QDict_ChlaProt_Path1_d3s7"

let headerLog = "StepN\tDirectionN\tGainBefore\tGainNow" 
let pathLog = @"c:\Users\mikha\Work-CSB\Projects\SSN\results\walkResults\"

File.AppendAllLines((sprintf "%s%s.txt" pathLog fileLogName), [headerLog])

let mutable stepCount = 0

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
                    let itemsChild = ids |> Array.map (fun i -> data.[i]) |> Array.concat |> groupIDFn
                    let itemsParent = data |> Array.concat |> groupIDFn 
                    getStepGainFn fn itemsChild itemsParent itemsParent.Length matrixSingles
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


type QDictionary_GValue =
    {
    NextStepGain: float
    NextStepState: int list list
    mutable MaxGain: float
    } 

let walkingFn kmeanKKZ depth matrixSingletons (singles: Map<string,Node<string,Item>>) gainFn (data: Map<string, Item []>) = 
    
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

        //let fileLogInit = sprintf "0.\t0.\tnan\tnan\t0.\t%f" (fst initialState)
        //File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogInit])

        let fileLog = sprintf "0.\t0.\t0.\t%f" (fst initialState)
        File.AppendAllLines((sprintf "%s%s.txt" pathLog fileLogName), [fileLog])

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
                && (countDirections<3) && (iStep<7) do // (pq.Top() > gainCurrent) && (countDirections<3) && (iStep<6)
                
                    countDirections <- countDirections + 1
                    
                    // order represents the moving: a - moved element, b - target cluster
                    let (a,b) = pairArray.[pq.TopIndex()]                 

                    //let fileLogStep = sprintf "%i\t%i\t%i\t%i\t%f\t%f" iStep countDirections a b gainCurrent (pq.Top())
                    //File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogStep])

                    let mA_new = snd mG.[a].[b]

                    let fileStep = sprintf "%i\t%i\t%f\t%f" iStep countDirections gainCurrent (pq.Top())
                    File.AppendAllLines((sprintf "%s%s.txt" pathLog fileLogName), [fileStep])

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

                    if (qDictionary |> Map.containsKey mA_new) then
                        
                        //let fileLogState = sprintf "%A was already visited, no step further" (reflectStateID mA_new)
                        //File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogState]) 
                        
                        yield! [] // how to get rid of the unnecessary empty lists? 
                    else
                    
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

                        stepCount <- stepCount + 1

                        yield (stats)
                        yield! loop (iStep+1) mA_new pairArrayNew matrixG pq_new new_moved
            ]
    
        Seq.appendSingleton (loop 1 matrixA pairArray' matrixG_origin pq_origin [||]) initialState // is it really neccessary to hand over pairArray?!
        

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
    

// incorporate the pathWalking into the createTree

type ClusterFn = int -> Map<string, Item []> -> (string * (Item [])) [] []

type WalkingFn = ClusterFn -> int -> float [,] -> Map<string,Node<string,Item>> -> (float -> float -> int -> int -> float) -> Map<string, Item []> -> Map<string, Item []>

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
                                                
                        //let fileLogData = sprintf "for elements:%A" (x |> Map.toArray |> Array.map fst)
                        //File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogData])

                        let fileLogData = sprintf "for elements:%A" (x |> Map.toArray |> Array.map fst)
                        File.AppendAllLines((sprintf "%s%s.txt" pathLog fileLogName), [fileLogData])
                        File.AppendAllLines((sprintf "%s%s.txt" pathLog fileLogName), [headerLog])

                        let singles = 
                            x
                            |> Map.fold (fun state key nodes ->
                                    let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrixSingletons
                                    state 
                                    |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)

                        walkFn kmeanKKZ depth matrixSingletons singles gainFn x 
                        |> Map.fold (fun state key nodes ->
                            match (singles.TryFind key) with
                            | None ->
                                let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn nodeMembers) matrixSingletons
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

//!!! Rethink idea of directions!!! not element -> element but element -> cluster

let applySSNcombi data setN = createTreeWalking (SSN.getStepGainNodeSetnR setN) None SSN_combi data
let applySSN data setN = createTreeWalking (SSN.getStepGainNodeSetnR setN) None (SSN_walk (FunctionsExp.PreClusterFunctions.kmeanGroupsKKZ, (walkingFn))) data
let readMM data setN = createTreeWalking (SSN.getStepGainNodeSetnR setN) None MM data

applySSN (TestData.SynData.data1) 10
applySSNcombi (TestData.SynData.data1) 10

stepCount

let dataPOI = 
    ChlamyProteome.dataAll
    |> Array.filter (fun x -> x.BinL.[0]="29")
    |> Array.mapi (fun id x -> {x with ID=id})

let data8 = 
    ArabiProteome.itemsWithMapManFound
    |> Array.filter (fun x -> x.BinL.[0]="8")
    |> Array.mapi (fun id x -> {x with ID=id})

let dataCut = [|data8.[0]; data8.[11]; data8.[75]; data8.[69]|] |> Array.mapi (fun id i -> {i with ID=id})

let treeMM = readMM data8 data8.Length
treeMM.GroupGain
let walk = applySSN dataPOI dataPOI.Length      
walk.GroupGain                              
let combi = applySSNcombi dataPOI dataPOI.Length
combi.GroupGain

let matrixOfN = 
    data8
    |> distMatrixWeightedOf distanceMatrixWeighted None
    /// normalised
let dissMax = 
    matrixOfN
    |> Array2D.array2D_to_seq
    |> Seq.max

let normMatrix = 
    matrixOfN
    |> Array2D.map (fun i -> i/dissMax)

walk
|> Tree.filterLeaves
|> Analysis.pointDxC normMatrix
|> fun (x,y) -> weightedEuclidean None [x;y] [0.;0.]

let timeAra = [|15.;180.;2880.;5760.;(5760.+15.);(5760.+180.);(5760.+2880.);(5760.+5760.)|]
Plots.drawKinetik (walk |> Tree.findNode ["1";"3"]) timeAra "subbin 8.3" |> Chart.Show//AsImage (StyleParam.ImageFormat.JPEG)

treeMM

GePhi.sendToGephiFromTreeParam treeMM

GePhi.sendToGephiFromTreeParam walk
GePhi.sendToGephiFromTreeParam combi


let pathsSorted =
    ArabiProteome.itemsWithMapManFound
    |> Array.groupBy (fun x -> x.BinL.[0])
    |> Array.filter (fun (bin,l) -> l.Length > 2)
    |> Array.sortBy (fun (bin,l) -> l.Length )
    |> Array.map (fst)
    
let pathsSortedChildren =
    ChlamyProteome.dataAll //ArabiProteome.itemsWithMapManFound
    |> Array.groupBy (fun x -> x.BinL.[0])
    |> Array.filter (fun (bin,l) -> l.Length > 2)
    |> Array.map (fun (bin,l) -> 
        let data = l |> Array.mapi (fun id x -> {x with ID=id})
        let tree = readMM data data.Length
        let childrenN = Tree.filterChildrenList tree
        (bin, childrenN |> List.max) )
    |> Array.filter (fun (bin,cn) -> cn<=10)
    |> Array.map fst

let pathFile = @"c:\Users\mikha\Work-CSB\Projects\SSN\results\walkResults\"
let neader = "Path\tcombi_GG\tcombi_Time\tcombi_DxC\twalk_GG\twalk_Time\twalk_DxC\tcomparison"

let lines = 
    [|"1";"31"|]
    |> List.ofArray
    |> List.map (fun path -> 
        let data = 
            ChlamyProteome.dataAll //ArabiProteome.itemsWithMapManFound
            |> Array.filter (fun x -> x.BinL.[0]=path)
            |> Array.mapi (fun id x -> {x with ID=id})
        let stopwatch = new System.Diagnostics.Stopwatch()
        stopwatch.Start()
        let walk = applySSN (data) data.Length
        let timeWalk = (stopwatch.Elapsed.TotalSeconds)
        let combi = applySSNcombi (data) data.Length
        let timeCombi = (stopwatch.Elapsed.TotalSeconds)
        stopwatch.Stop()
        let matrixOfN = 
            data
            |> distMatrixWeightedOf distanceMatrixWeighted None
            /// normalised
        let dissMax = 
            matrixOfN
            |> Array2D.array2D_to_seq
            |> Seq.max

        let normMatrix = 
            matrixOfN
            |> Array2D.map (fun i -> i/dissMax)

        let dxcW = 
            walk
            |> Tree.filterLeaves
            |> Analysis.pointDxC normMatrix
            |> fun (x,y) -> weightedEuclidean None [x;y] [0.;0.]
    
        let dxcC = 
            combi
            |> Tree.filterLeaves
            |> Analysis.pointDxC normMatrix
            |> fun (x,y) -> weightedEuclidean None [x;y] [0.;0.]

        printfn "path: %s, %i" path data.Length 
        printfn "combi GG; combi Time; walk GG; walk Time; walk DxC; comparison"
        printfn "%s %f %f %f %f %f %f %i" path combi.GroupGain timeCombi dxcC walk.GroupGain timeWalk dxcW (Tree.treeComparison combi walk)
        sprintf "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%i" path combi.GroupGain timeCombi dxcC walk.GroupGain timeWalk dxcW (Tree.treeComparison combi walk)

        //(path, data.Length, combi.GroupGain, timeCombi, dxcC, walk.GroupGain, timeWalk, dxcW, (Tree.treeComparison combi walk))
        )


let linesWAcombi = 
    [|"29"|]
    |> List.ofArray
    |> List.map (fun path -> 
        let data = 
            ChlamyProteome.dataAll //ArabiProteome.itemsWithMapManFound
            |> Array.filter (fun x -> x.BinL.[0]=path)
            |> Array.mapi (fun id x -> {x with ID=id})
        let stopwatch = new System.Diagnostics.Stopwatch()
        stopwatch.Start()
        let walk = applySSN (data) data.Length
        let timeWalk = (stopwatch.Elapsed.TotalSeconds)
        //let combi = applySSNcombi (data) data.Length
        //let timeCombi = (stopwatch.Elapsed.TotalSeconds)
        stopwatch.Stop()
        let matrixOfN = 
            data
            |> distMatrixWeightedOf distanceMatrixWeighted None
            /// normalised
        let dissMax = 
            matrixOfN
            |> Array2D.array2D_to_seq
            |> Seq.max

        let normMatrix = 
            matrixOfN
            |> Array2D.map (fun i -> i/dissMax)

        let dxcW = 
            walk
            |> Tree.filterLeaves
            |> Analysis.pointDxC normMatrix
            |> fun (x,y) -> weightedEuclidean None [x;y] [0.;0.]
    
        //let dxcC = 
        //    combi
        //    |> Tree.filterLeaves
        //    |> Analysis.pointDxC normMatrix
        //    |> fun (x,y) -> weightedEuclidean None [x;y] [0.;0.]

        printfn "path: %s, %i" path data.Length 
        printfn "combi GG; combi Time; walk GG; walk Time; walk DxC; comparison"
        printfn "%s %f %f %f %f %f %f %i" path 0. 0. 0. walk.GroupGain timeWalk dxcW 0
        sprintf "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%i" path 0. 0. 0. walk.GroupGain timeWalk dxcW 0

        //(path, data.Length, combi.GroupGain, timeCombi, dxcC, walk.GroupGain, timeWalk, dxcW, (Tree.treeComparison combi walk))
        )

File.AppendAllLines((sprintf "%s%s.txt" pathFile "oldChlamyProtData_"), neader :: linesWAcombi)


/////////////////////////// ####### plot path walking tree

open FSharpAux.IO

type DoubleArrayConverter() = 
        inherit SchemaReader.Attribute.ConverterAttribute()
        override this.convertToObj = 
            SchemaReader.Converter.Collection (fun (strs : seq<string>) -> 
                                                    (strs |> Seq.map (fun s -> FSharpAux.String.tryParseFloatDefault nan s) |> Seq.toArray) |> box)

type ProteinItemRead = {
        [<SchemaReader.Attribute.FieldAttribute([| "StepN";"DirectionN";"movedID";"MovedTo";"GainBefore";"GainNow" |])>]  [<DoubleArrayConverter>]    Features        : float []
        }

//// Read Proteins database
let reader = new SchemaReader.Csv.CsvReader<ProteinItemRead>(schemaMode=SchemaReader.Csv.SchemaMode.Fill)

/// Variable, contains all raw data from protein DataBase file in csvPath
let dataProtein file = 
    reader.ReadFile (file, General.separatorTab, General.hasHeader) 
    |> Seq.toArray

let linesFromFile = 
    dataProtein @"C:\Users\mikha\Work-CSB\Projects\SSN\results\walkResults\path8-2-11_k2_toPlot_problematic.txt"
    |> Array.map (fun x -> [(x.Features.[0]-1.,x.Features.[4]);(x.Features.[0],x.Features.[5])],x.Features.[1])

linesFromFile.Length

let closestPath (lines: (((float * float) list ) * float) [] ) = 
    let max = lines |> Array.maxBy (fun (i,_) -> snd i.[1]) |> fst |> List.item 1 |> snd
    let closOpt = lines |> Array.filter (fun (i,_) -> (snd i.[1]) = max ) |> Array.minBy (fun (i,_) -> fst i.[1])

    let rec loop (nextStep: ((float * float) list ) * float) = 
        [
        let prevStep = lines |> Array.find (fun (i,_) -> (fst i.[1])=(fst (fst nextStep).[0]) && (snd i.[1])=(snd (fst nextStep).[0]) )
        if ((fst prevStep).[0] |> fst) = 0. then
            yield prevStep
        else
            yield prevStep
            yield! loop prevStep]

    closOpt :: (loop closOpt)

let plotGainWalk (a: (((float * float) list ) * float) [] ) =
    let color = function
        |1 -> colorBlue
        |2 -> colorBrightBlue
        |3 -> colorGreen
        |4 -> colorOrange
        |5 -> colorYellow
        |_ -> colorGray
 
    let lines = 
        a
        |> Array.map (fun (i,j) -> 
                        let c = color (int j)
                        Chart.Line (i, Color = c)
                    )
        |> Array.toList
        |> List.rev

    let pathToOptimal = 
        a
        |> closestPath
        |> (fun l ->
            (l |> List.head |> fst).[1]::(l |> List.map (fun (i,_) -> i.[0]))
            )
                        
    (Chart.Line (pathToOptimal, Color = "rgba(0,0,0,1)"))::lines |> List.rev |> Chart.Combine |> Chart.Show

plotGainWalk linesFromFile



closestPath linesFromFile
|> (fun l ->
            (l |> List.head |> fst).[1]::(l |> List.map (fun (i,_) -> i.[0]))
            )