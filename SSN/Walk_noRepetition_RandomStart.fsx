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

open System 
open System.IO 
open FSharpAux
open FSharp.Plotly
open FSharp.Stats

open BioFSharp

open Functions
open Functions.General
open Functions.Walk
open TestData
open GePhi
open Types
open PQ
open Auxilliary
open Plots
//open FunctionsExp

#time

let mutable stepCount = 0

let delta = 0.06

let walkingFn kmeanKKZ depth matrixSingletons (singles: Map<string,Node<string,Item>>) gainFn (dataM: Map<string, Item []>) = 
    
    let mutable qDictionary: Map<(int list list),((float*(int list list)) [] [])> = Map.empty

    let dataGroupsA = dataM |> Map.toArray

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
                && (pq.Top()>(gainCurrent - delta*gainCurrent)) // no sinking lower delta
                && (countDirections < 1) // max direction checked = 1
                && (iStep < (dataM.Count/1)) //&& (iStep<5) // max path length = 5
                do 
                
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

                        stepCount <- stepCount + 1

                        yield (stats)
                        yield! loop (iStep+1) mA_new pairArrayNew matrixG pq_new new_moved
            ]
    
        Seq.appendSingleton (loop 1 matrixA pairArray' matrixG_origin pq_origin [||]) initialState 

    (seq [singleGG |> Array.sum, Array.init dataM.Count (fun i -> [|i|])]) :: 
        ([2 .. (dataM.Count-1)] 
        |> List.map (fun i -> 

            //let fileLogKMEAN = sprintf "Kmean with k=%i"  i
            //File.AppendAllLines((sprintf "%s%s.txt" pathLOG fileSubName), [fileLogKMEAN])

            dataM
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

let walkingHC_Fn hcFn depth matrixSingletons (singles: Map<string,Node<string,Item>>) gainFn (dataM: Map<string, Item []>) = 
    
    let mutable qDictionary: Map<(int list list),((float*(int list list)) [] [])> = Map.empty

    let dataGroupsA = dataM |> Map.toArray

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
                (pq.Length > 0) 
                && (pq.Top() > (gainCurrent - delta*gainCurrent) ) // no sinking lower delta
                && (countDirections < 1) // max direction checked = 1
                && (iStep < (dataM.Count/1) ) //&& (iStep<5) // max path length = 5
                do 
                
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

                        stepCount <- stepCount + 1

                        yield (stats)
                        yield! loop (iStep+1) mA_new pairArrayNew matrixG pq_new new_moved
            ]
    
        Seq.appendSingleton (loop 1 matrixA pairArray' matrixG_origin pq_origin [||]) initialState 

    //let paralFn i =
    //    async {return (i |> superFunctionTestG dataGroups singleGG gainFn matrixSingletons)}

    //(seq [singleGG |> Array.sum, Array.init dataM.Count (fun i -> [|i|])]) :: 
    //    ([2 .. 4 .. (dataM.Count-1)] 
    //    |> hcFn dataM
    //    |> List.toArray
    //    |> Array.map (fun i -> paralFn i)
    //    |> Async.Parallel
    //    |> Async.RunSynchronously
    //    |> List.ofArray
    //    )

    let eval fn matrixSingles (initConf: (string * Item [] ) [] [])  =
        (initConf
        |> Array.sumBy (fun cluster ->
            if cluster.Length=1 then 
                (Map.find (fst cluster.[0]) singles).GroupGain
            else
                let itemsChild = cluster |> Array.map (snd) |> Array.concat |> General.groupIDFn
                let itemsParent = dataGroups |> Array.concat |> General.groupIDFn 
                General.getStepGainFn fn itemsChild itemsParent itemsParent.Length matrixSingles
            ), initConf)

    (seq [singleGG |> Array.sum, Array.init dataM.Count (fun i -> [|i|])]) :: 
        ([2 .. (dataM.Count-1)] 
        |> hcFn dataM
        |> List.sortByDescending ((eval gainFn matrixSingletons) >> fst)
        |> fun x -> if x.Length>10 then x.[0 .. 10] else x
        |> List.map (fun i -> i |> superFunctionTestG dataGroups singleGG gainFn matrixSingletons)
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
    
let walkingParalFn kmeanKKZ depth matrixSingletons (singles: Map<string,Node<string,Item>>) gainFn (dataM: Map<string, Item []>) = 
    
    let mutable qDictionary: Map<(int list list),((float*(int list list)) [] [])> = Map.empty

    let dataGroupsA = dataM |> Map.toArray

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
                (pq.Length > 0) 
                && (pq.Top() > gainCurrent)// - delta*gainCurrent)) // no sinking lower delta
                && (countDirections < 1) // max direction checked = 1
                && (iStep < 5)//(dataM.Count/4)) //&& (iStep<5) // max path length = 5
                do 
                
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

                        stepCount <- stepCount + 1

                        yield (stats)
                        yield! loop (iStep+1) mA_new pairArrayNew matrixG pq_new new_moved
            ]
    
        Seq.appendSingleton (loop 1 matrixA pairArray' matrixG_origin pq_origin [||]) initialState 

    let paralFn i =
        async {return (dataM |> kmeanKKZ i |> superFunctionTestG dataGroups singleGG gainFn matrixSingletons)}

    (seq [singleGG |> Array.sum, Array.init dataM.Count (fun i -> [|i|])]) :: 
        ([|2 .. 2 .. (dataM.Count-1)|] 
        |> Array.map (Array.create 10)
        |> Array.concat
        |> Array.map (fun i -> paralFn i)
        |> Async.Parallel
        |> Async.RunSynchronously
        |> List.ofArray
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

/// set a random start
let randomConf (items: Map<string,Item []>) =
    let n = items.Count
    let r = new System.Random()
    
    let k = r.Next(2,n)

    items
    |> Map.toArray
    |> Array.shuffleFisherYates  
    |> Array.indexed
    |> Array.groupBy (fun (i,_) -> if i<k then i else r.Next(0,k) )
    |> Array.map (snd >> Array.map (snd))
    
//randomConf (breakGroup TestData.SynData.data1 0)

/// set a random start with given k clusters
let randomKConf k (items: Map<string,Item []>) =

    let r = new System.Random()
    
    items
    |> Map.toArray
    |> Array.shuffleFisherYates  
    |> Array.indexed
    |> Array.groupBy (fun (i,_) -> if i<k then i else r.Next(0,k) )
    |> Array.map (snd >> Array.map (snd))
    
//randomKConf 3 (breakGroup TestData.SynData.data1 0)

//// apply

let applySST_clustDet setN data = applySST_onlyHierClust setN data

let applySST_clust setN data = applySST_onlyKMClust setN data


let applySST_walk_write setN data = 
    (createTree (getStepGainNodeSetnR setN) None (SST_walk (Walk.walkingFn 1 5 Clustering.kmeanGroupsKKZ)) data)
        
let applySST_walkRandom_write setN data = 
    (createTree (getStepGainNodeSetnR setN) None (SST_walk (walkingFn randomKConf)) data)

let applySST_walkFastRandom setN data = 
    (createTree (getStepGainNodeSetnR setN) None (SST_walk (walkingParalFn randomKConf)) data)

let applySST_walkFromHC setN data = 
    (createTree (getStepGainNodeSetnR setN) None (SST_walk  (walkingHC_Fn Clustering.clusterHierGroups)) data)

let asyncApplySSN data setN = async {return (applySST_walk setN data )}


let dataPOI = 
    ArabiProteome.itemsWithMapManFound
    |> Array.filter (fun x -> x.BinL.[0]="8")
    |> Array.mapi (fun id x -> {x with ID=id})
dataPOI.Length

let tree1 = applySST_walkRandom_write (dataPOI.Length) dataPOI
let tree1_ = applySST_walkFastRandom (dataPOI.Length) dataPOI

let tree1_hc_walk_ViertelK = applySST_walkFromHC (dataPOI.Length) dataPOI

let tree1_hc = applySST_clustDet (dataPOI.Length) dataPOI 

let tree8_hc_walk = applySST_walkFromHC (dataPOI.Length) dataPOI


tree1.GroupGain         // 69.4115
tree1_.GroupGain        // 48.6404
tree1_hc.GroupGain      // 71.6172 good enough?
tree8_hc_walk.GroupGain // 79.2476 for 1.5 hr
//tree1_hc_walk_halfK.GroupGain // 79.2476 for 3 min!!!
tree1_hc_walk_ViertelK.GroupGain // 79.2476 for 40 sec!!! or 1 min with parallelization.... weird.

//sendToGephiFromTreeParam tree1_hc_walk

let subbin_1_3 = (tree1 |> Tree.findNode ["1";"3"]) |> Tree.filterLeaves

Plots.drawKinetikRange [|0 .. 7|] "path 1.3 for ArabiProteome" subbin_1_3 |> Chart.Show

let dataSet =
    [|ChlamyProteome.dataAll
    |> Array.filter (fun x -> x.BinL.[0]="1")
    |> Array.mapi (fun id x -> {x with ID=id});
    ChlamyProteome.dataAll
    |> Array.filter (fun x -> x.BinL.[0]="13")
    |> Array.mapi (fun id x -> {x with ID=id});
    ChlamyProteome.dataAll
    |> Array.filter (fun x -> x.BinL.[0]="11")
    |> Array.mapi (fun id x -> {x with ID=id});
    ChlamyProteome.dataAll
    |> Array.filter (fun x -> x.BinL.[0]="31")
    |> Array.mapi (fun id x -> {x with ID=id})|]

let trees =
    dataSet
    |> Array.map (fun x -> asyncApplySSN x x.Length)
    |> Async.Parallel
    |> Async.RunSynchronously

trees.[3]

let pathsSorted =
    ArabiProteome.itemsWithMapManFound
    |> Array.groupBy (fun x -> x.BinL.[0])
    |> Array.filter (fun (bin,l) -> l.Length > 2)
    |> Array.map (fun (bin,l) -> 
        let data = l |> Array.mapi (fun id x -> {x with ID=id})
        let tree = readMM data.Length data
        let childrenN = Tree.filterChildrenList tree
        printfn "%s\t%i\t%i\t%f" bin data.Length (childrenN |> List.max) tree.GroupGain
        (bin, data.Length, childrenN |> List.max, tree.GroupGain) )
    
let pathsSortedChildren =
    ArabiProteome.itemsWithMapManFound  // ChlamyProteome.dataAll //   ArabiTranscriptome.itemsWithMapManIdentified //
    |> Array.groupBy (fun x -> x.BinL.[0])
    |> Array.filter (fun (bin,l) -> l.Length > 2 && bin<>"35")
    |> Array.map (fun (bin,l) -> 
        let data = l |> Array.mapi (fun id x -> {x with ID=id})
        let tree = readMM data.Length data
        let childrenN = Tree.filterChildrenList tree
        (bin, childrenN |> List.max, tree.GroupGain) )
    //|> Array.filter (fun (_,cn,_) -> cn<11)
    |> Array.sortBy (fun (_,cn,_) -> cn)
    |> Array.map (fun (x,_,_) -> x)


let pathFile = sprintf @"%sresults\walk_hcStart\" General.pathToData
let neader = "Path\tnRoot\tchildrenMax\thc_GG\thc_Time\thc_DxC\twalk_GG\twalk_Time\twalk_DxC"

File.AppendAllLines((sprintf "%s%s.txt" pathFile "ArabiProtData_hc10_d1_sSize_del6"), [neader])
       
let linesWAcombi =
    pathsSortedChildren // [|"22"; "6"; "24"; "5"; "14"; "17"; "2"; "25"; "8"; "12"; "7"; "23"; "16"; "18"; "15"; "13"; "10"|] //[|"13"; "9"; "11"; "4"; "34"; "3"; "19"; "28"; "21"; "26"; "33"; "1"; "27"; "20"; "30"; "29"; "31";|]
    |> List.ofArray
    |> List.map (fun path ->
        let data =
            ArabiProteome.itemsWithMapManFound //  ChlamyProteome.dataAll //  
            |> Array.filter (fun x -> x.BinL.[0]=path)
            |> Array.mapi (fun id x -> {x with ID=id})
        let tree = readMM  data.Length data
        let childrenN = Tree.filterChildrenList tree |> List.max
        let stopwatch = new System.Diagnostics.Stopwatch()
        stopwatch.Start()
        let hc = applySST_clustDet data.Length data // applySSNcombi data.Length data // 
        let timeHC = (stopwatch.Elapsed.TotalSeconds)
        stopwatch.Stop()
        stopwatch.Start()
        let walk = applySST_walkFromHC data.Length data
        let timeWalk = (stopwatch.Elapsed.TotalSeconds)
        stopwatch.Stop()

        let matrixOfN =
            data
            |> General.distMatrixWeightedOf General.distanceMatrixWeighted None
           
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
            |> fun (x,y) -> General.weightedEuclidean None [x;y] [0.;0.]

        let dxcHC =
            hc
            |> Tree.filterLeaves
            |> Analysis.pointDxC normMatrix
            |> fun (x,y) -> General.weightedEuclidean None [x;y] [0.;0.]

        printfn "path; nRoot; childrenMax; hc GG; hc Time; hc DxC; walk GG; walk Time; walk DxC"
        printfn "%s %i %i %f %f %f %f %f %f"
            path data.Length childrenN hc.GroupGain timeHC dxcHC walk.GroupGain timeWalk dxcW 

        let treeLinesW =
            walk
            |> Tree.filterLeaves
            |> Array.concat
            |> Array.map (fun x ->  (sprintf "%i" x.ID) + "\t" + (String.Join(";", x.OriginalBin))  + "\t" + (String.Join(";", x.BinL)))
        File.AppendAllLines(sprintf @"c:\Users\mikha\source\repos\SSN\results\walk_hcStart\Arabi_path_walk_hc10best_%s.txt" path, treeLinesW)

        //let treeLinesHC =
        //    hc
        //    |> Tree.filterLeaves
        //    |> Array.concat
        //    |> Array.map (fun x ->  (sprintf "%i" x.ID) + "\t" + (String.Join(";", x.OriginalBin))  + "\t" + (String.Join(";", x.BinL)))
        //File.AppendAllLines(sprintf @"c:\Users\mikha\source\repos\SSN\results\walk_hcStart\Chlamy_hc_%s.txt" path, treeLinesHC)

        let line =
            sprintf "%s\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f"
                path data.Length childrenN hc.GroupGain timeHC dxcHC walk.GroupGain timeWalk dxcW 
       
        File.AppendAllLines((sprintf "%s%s.txt" pathFile "ArabiProtData_hc10_d1_sSize_del6"), [line])

        line
        )

        
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
    dataProtein @"..\results\path8-2-11_k2_toPlot_problematic.txt"
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