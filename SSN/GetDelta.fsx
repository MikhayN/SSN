#r @"..\lib\MathNet.Numerics.dll"
#r @".\bin\Debug\netstandard.dll"

#r @".\bin\Debug\BioFSharp.dll"
#r @".\bin\Debug\BioFSharp.IO.dll"
#r @".\bin\Debug\FSharpAux.dll"
#r @".\bin\Debug\FSharpAux.IO.dll"
#r @".\bin\Debug\Newtonsoft.Json.dll"
#r @".\bin\Debug\FSharp.Stats.dll"
#r @".\bin\Debug\FSharp.Plotly.dll"
#r @".\bin\Debug\FSharpGephiStreamer.dll"

#load "Types.fs"
#load "PQ.fs"
#load "Functions.fs"
#load "FunctionsExp.fs"
#load "GePhi.fs"
#load "TestData.fs"
#load "SOM.fs"
#load "Auxilliary.fs"
#load "Plots.fs"
#load "ElbowCriteria.fs"

open System 
open System.IO 
open FSharpAux
open FSharpAux.IO
open FSharp.Plotly
open FSharp.Stats

//open Functions
open Functions
open Functions.SSN
open Functions.Walk
open TestData
open GePhi
open Types
open PQ
open Auxilliary
open Plots
open ElbowCriteria

#time

///////////////////////////

//let fileLogName = "QDict_ChlaProt_Path29_d2s5"

let headerLog = "StepN\tDirN\tGainBefore\tGainNow" 
let pathLog = sprintf @"%sresults\Arabi\" General.pathToData

//File.AppendAllLines((sprintf "%s%s.txt" pathLog fileLogName), [headerLog])

let mutable countOverall = 0

let walkingFnWrite fileLogName kmeanKKZ depth matrixSingletons (singles: Map<string,Node<string,Item>>) gainFn (data: Map<string, Item []>)  = 
    
    //let fileNameSubindex = sprintf "%i-%i" depth ((data |> Map.toArray).[0] |> snd).[0].ID
    let fileNameSubindex = sprintf "%i" countOverall
    File.AppendAllLines((sprintf "%s%s\\%s.txt" pathLog fileLogName fileNameSubindex), [headerLog])
    //File.AppendAllLines((sprintf "%sQDict_%s.txt" pathLog fileLogName), ["new walking started"; headerLog])

    countOverall <- countOverall + 1

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

        let initStateIDs = reflectStateID matrixA
        let initialState = (gainFnn data singleGG matrixSingles fn initStateIDs, initStateIDs)

        let fileLog = sprintf "0\t0\t0.\t%f" (fst initialState)
        //File.AppendAllLines((sprintf "%sQDict_%s_%s.txt" pathLog fileLogName fileNameSubindex), [fileLog])
        File.AppendAllLines((sprintf "%s%s\\%s.txt" pathLog fileLogName fileNameSubindex), [fileLog])

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
                && (countDirections<1) && (iStep<5) do // (pq.Top() > gainCurrent) && (countDirections<3) && (iStep<6)
                
                    countDirections <- countDirections + 1
                    
                    // order represents the moving: a - moved element, b - target cluster
                                     
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

                        let fileStep = sprintf "%i\t%i\t%f\t%f" iStep countDirections gainCurrent (pq.Top())
                        //File.AppendAllLines((sprintf "%sQDict_%s_%s.txt" pathLog fileLogName fileNameSubindex), [fileStep])
                        File.AppendAllLines((sprintf "%s%s\\%s.txt" pathLog fileLogName fileNameSubindex), [fileStep])
                        
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

///////////////////////////

let applySST_walk_write setN path data = 
    countOverall <- 0
    (SSN.createTree (SSN.getStepGainNodeSetnR setN) None (SST_walk (Clustering.kmeanGroupsKKZ, (walkingFnWrite path) )) data)

let dataSet =
    ArabiProteome.itemsWithMapManIdentified //ChlamyProteome.dataAll
    |> Array.groupBy (fun x -> x.BinL.[0])
    |> Array.map (snd >> Array.mapi (fun id x -> {x with ID=id}))
    |> Array.filter (fun il -> il.Length>2 && il.[0].OriginalBin.[0]<>"35")
    |> Array.sortBy (fun data -> 
        let tree = readMM data.Length data 
        let childrenN = Tree.filterChildrenList tree
        childrenN |> List.max
    )

dataSet |> Array.iter (fun i -> printfn "%s" i.[0].OriginalBin.[0])


let dataPOI = 
    ArabiProteome.itemsWithMapManIdentified //ChlamyProteome.dataAll
    |> Array.filter (fun x -> x.BinL.[0]="23")
    |> Array.sortBy (fun x -> x.ProteinL)
    |> Array.mapi (fun id x -> {x with ID=id})

applySST_walk_write dataPOI.Length (dataPOI.[0].OriginalBin.[0]) dataPOI


let treesTracking =
    dataSet
    |> Array.map (fun x -> 
        printfn "starting bin %s" x.[0].OriginalBin.[0]
        applySST_walk_write x.Length (x.[0].OriginalBin.[0]) x)



/////////////////////////// ####### plot path walking tree

type ProteinItemRead = {
        [<SchemaReader.Attribute.FieldAttribute("StepN" )>]     [<ArabiTranscriptome.DoubleConverter>]    StepN        : float 
        [<SchemaReader.Attribute.FieldAttribute("DirN")>]       [<ArabiTranscriptome.DoubleConverter>]    DirN         : float 
        [<SchemaReader.Attribute.FieldAttribute("GainBefore")>] [<ArabiTranscriptome.DoubleConverter>]    GainBefore   : float 
        [<SchemaReader.Attribute.FieldAttribute("GainNow" )>]   [<ArabiTranscriptome.DoubleConverter>]    GainNow      : float 
        }

//// Read Proteins database
let reader = new SchemaReader.Csv.CsvReader<ProteinItemRead>(schemaMode=SchemaReader.Csv.SchemaMode.Fill)

/// Variable, contains all raw data from protein DataBase file in csvPath
let dataProtein file = 
    reader.ReadFile (file, General.separatorTab, true) 
    |> Seq.toArray

let linesFromFile file = 
    dataProtein (sprintf @"%s\results\%s.txt" General.pathToData file)
    //|> Array.map (fun x -> [(x.Features.[0]-1.,x.Features.[2]);(x.Features.[0],x.Features.[3])],x.Features.[1])

let tryOn = (linesFromFile "Arabi\\8\\1")//.Length

let closestPath (lines: ProteinItemRead [] ) = 

    let closOpt = (lines |> Array.maxBy (fun i -> i.GainNow))
    
    let rec loop (nextStep: ProteinItemRead) = 
        [
        if nextStep.GainBefore = 0. then
            yield nextStep
        else
            yield nextStep
            let prevStep = lines |> Array.find (fun (i) -> i.GainNow = nextStep.GainBefore )
            yield! loop prevStep]

    (loop closOpt)

closestPath tryOn

//////////////////

let startMinEnd (closestPath: ProteinItemRead list ) =
    let endG = (closestPath.Head).GainNow
    let start =  (closestPath |> List.last).GainNow
    let maxDeltaG = closestPath |> List.maxBy (fun x -> x.GainBefore - x.GainNow ) // if negative, no sinking
    (start, maxDeltaG, endG)

startMinEnd (closestPath tryOn)

let getSingleDelta (sme: float * ProteinItemRead * float) =
    let (s,mD,e) = sme
    (mD.GainBefore - mD.GainNow)/(mD.GainBefore) // if negative, don't count

let averageDelta (dataD: (float*ProteinItemRead*float) list) =
    if dataD=[] then 0.
    else
        let listChosen = 
            dataD
            |> List.choose 
                (fun (s,mD,e) -> 
                    if e>=s then
                        let delta = getSingleDelta (s,mD,e)
                        if delta>0. then
                            Some delta
                        else None
                    else None)
        if listChosen = [] then 0.
        else
            listChosen
            |> List.average


let maxDelta (dataD: (float*ProteinItemRead*float) list) =
    if dataD=[] then 0.
    else
        let listChosen = 
            dataD
            |> List.choose 
                (fun (s,mD,e) -> 
                    if e>=s then
                        let delta = getSingleDelta (s,mD,e)
                        if delta>0. then
                            Some delta
                        else None
                    else None)
        if listChosen = [] then 0.
        else
            listChosen
            |> List.max

let manyPaths bin maxFileName =
    [0 .. maxFileName]
    |> List.map (fun i -> 
        i
        |> sprintf "Arabi\\%i\\%i" bin
        |> linesFromFile)
    |> List.filter (fun i -> i<>[||])
    |> List.map
        (closestPath >> startMinEnd)
    |> maxDelta

// {number of step, avoided by setting delta} against {number of steps, taken despite the delta}
let estimLowerThanDelta bin maxFileName delta =
    let parted =
        [|0 .. maxFileName|]
        |> Array.map (fun i -> 
            i
            |> sprintf "Arabi\\%i\\%i" bin
            |> linesFromFile)
        |> Array.filter (fun i -> i<>[||])
        |> Array.concat
        |> Array.partition (fun i -> (i.GainBefore - i.GainNow)/i.GainBefore >= delta)
    //parted
    ( float (fst parted).Length ) / ( float (snd parted).Length )

// Chlamy
let deltas29 = manyPaths 29 41

// Arabi -> set { DELTA = 6% }
let deltas8 = manyPaths 8 19 // 3% aver // 5.7% max
let deltas11 = manyPaths 11 25 // no sinking 
let deltas23 = manyPaths 23 21 // no sinking
let deltas22 = manyPaths 22 1 // no sinking
let deltas6 = manyPaths 6 2 // no sinking
let deltas24 = manyPaths 24 2 // no sinking
let deltas5 = manyPaths 5 2 // no sinking
let deltas14 = manyPaths 14 3 // no sinking
let deltas17 = manyPaths 17 22 // 0.3% aver // 0.3% max
let deltas2 = manyPaths 2 22 // no sinking
let deltas25 = manyPaths 25 5 // no sinking
let deltas12 = manyPaths 12 6 // no sinking
let deltas7 = manyPaths 7 4 // no sinking
let deltas16 = manyPaths 16 37 // no sinking
let deltas18 = manyPaths 18 7 // no sinking
let deltas15 = manyPaths 15 2 // no sinking
let deltas13 = manyPaths 13 64 // 4.6% aver // 4.6% max

let estim8 = estimLowerThanDelta 8 19 0.06 // 37 vs 100
let estim11 = estimLowerThanDelta 11 25 0.06 // 24 vs 100
let estim23 = estimLowerThanDelta 23 21 0.06 // 39 vs 100
let estim22 = estimLowerThanDelta 22 1 0.06 // 50 vs 100
let estim6 = estimLowerThanDelta 6 2 0.06 // 60 vs 100
let estim24 = estimLowerThanDelta 24 2 0.06 // 15 vs 100
let estim5 = estimLowerThanDelta 5 2 0.06 // 19 vs 100
let estim14 = estimLowerThanDelta 14 3 0.06 // 52 vs 100
let estim17 = estimLowerThanDelta 17 22 0.06 // 40 vs 100
let estim2 = estimLowerThanDelta 2 22 0.06 // 36 vs 100
let estim25 = estimLowerThanDelta 25 5 0.06 // 22 vs 100
let estim12 = estimLowerThanDelta 12 6 0.06 // 10 vs 100
let estim7 = estimLowerThanDelta 7 4 0.06 // 17 vs 100
let estim16 = estimLowerThanDelta 16 37 0.06 // 27 vs 100
let estim18 = estimLowerThanDelta 18 7 0.06 // 5 vs 100
let estim15 = estimLowerThanDelta 15 2 0.06 // 5 vs 100
let estim13 = estimLowerThanDelta 13 64 0.06 // 40 vs 100

/////////

let plotGainWalk (a: ProteinItemRead [] ) =
    let color = function
        |1 -> colorBlue
        |2 -> colorBrightBlue
        |3 -> colorGreen
        |4 -> colorOrange
        |5 -> colorYellow
        |_ -> colorGray
 
    let lines = 
        a
        |> Array.map (fun (i) -> 
                        let c = color (int i.DirN)
                        let line = [(i.StepN-1.,i.GainBefore);(i.StepN,i.GainNow)]
                        Chart.Line (line, Color = c)
                    )
        |> Array.toList
        |> List.rev

    let pathToOptimal = 
        a
        |> closestPath
        |> (fun l ->
            ((l |> List.head).StepN,(l |> List.head).GainNow) :: (l |> List.map (fun (i) -> (i.StepN-1.,i.GainBefore)))
            )
                        
    (Chart.Line (pathToOptimal, Color = "rgba(0,0,0,1)"))::lines |> List.rev |> Chart.Combine |> Chart.Show

plotGainWalk (tryOn)
