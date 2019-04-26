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

let headerLog = "StepN\tDirectionN\tGainBefore\tGainNow" 
let pathLog = sprintf @"%sresults\" General.pathToData

//File.AppendAllLines((sprintf "%s%s.txt" pathLog fileLogName), [headerLog])

let walkingFnWrite fileLogName kmeanKKZ depth matrixSingletons (singles: Map<string,Node<string,Item>>) gainFn (data: Map<string, Item []>)  = 
    
    let fileNameSubindex = sprintf "%i-%i" depth ((data |> Map.toArray).[0] |> snd).[0].ID
    //File.AppendAllLines((sprintf "%sQDict_%s.txt" pathLog fileLogName), ["new walking started"; headerLog])

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
        File.AppendAllLines((sprintf "%sQDict_%s_%s.txt" pathLog fileLogName fileNameSubindex), [fileLog])

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
                        File.AppendAllLines((sprintf "%sQDict_%s_%s.txt" pathLog fileLogName fileNameSubindex), [fileStep])

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

let applySST_walk_write setN path data = async {
    return (SSN.createTree (SSN.getStepGainNodeSetnR setN) None (SST_walk (Clustering.kmeanGroupsKKZ, (walkingFnWrite path) )) data)
    }

let dataSet =
    ChlamyProteome.dataAll
    |> Array.groupBy (fun x -> x.BinL.[0])
    |> Array.map (snd >> Array.mapi (fun id x -> {x with ID=id}))
    |> Array.filter (fun il -> il.Length>2 && il.[0].OriginalBin.[0]<>"35")

dataSet |> Array.iter (fun i -> printfn "%s" i.[0].OriginalBin.[0])


let treesTracking =
    dataSet
    |> Array.map (fun x -> applySST_walk_write x.Length (x.[0].OriginalBin.[0]) x)
    |> Async.Parallel
    |> Async.RunSynchronously


/////////////////////////// ####### plot path walking tree

open FSharpAux.IO

type DoubleArrayConverter() = 
        inherit SchemaReader.Attribute.ConverterAttribute()
        override this.convertToObj = 
            SchemaReader.Converter.Collection (fun (strs : seq<string>) -> 
                                                    (strs |> Seq.map (fun s -> FSharpAux.String.tryParseFloatDefault nan s) |> Seq.toArray) |> box)

type ProteinItemRead = {
        [<SchemaReader.Attribute.FieldAttribute([| "StepN";"DirectionN";"GainBefore";"GainNow" |])>]  [<DoubleArrayConverter>]    Features        : float []
        }

//// Read Proteins database
let reader = new SchemaReader.Csv.CsvReader<ProteinItemRead>(schemaMode=SchemaReader.Csv.SchemaMode.Fill)

/// Variable, contains all raw data from protein DataBase file in csvPath
let dataProtein file = 
    reader.ReadFile (file, General.separatorTab, General.hasHeader) 
    |> Seq.toArray

let linesFromFile path = 
    dataProtein (sprintf @"%s\results\WalkingSteps for Chlamy Proteome\QDict_%s.txt" General.pathToData path)
    |> Array.map (fun x -> [(x.Features.[0]-1.,x.Features.[2]);(x.Features.[0],x.Features.[3])],x.Features.[1])

(linesFromFile "28").Length

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

let startMinEnd (closestPath: ((float * float) list * float) list ) =
    let start = (closestPath.Head |> fst ).[0] |> snd
    let endG =  (closestPath |> List.last |> fst ).[1] |> snd
    let minG = closestPath |> List.map (fun x -> ((x |> fst).[1] |> snd) - ((x |> fst).[0] |> snd) ) |> List.min // if negative, no sinking
    (start, minG, endG)

let getSingleDelta (sme: float * float * float) =
    let (s,m,e) = sme
    (e-s)/(m) // if negative, don't count

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