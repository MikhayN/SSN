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
open FSharpAux.IO

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

//
////--------------------------------------------------------------------- Read Hermit classes ------------------------------------------------------------------
//

type JGIMapping =
    {
    JGI43full           : string
    JGI55full           : string
    JGI43               : string
    JGI55               : string
    MapManGMM           : string []
    MapManDescription   : string []
    }
    
let createJGIMapping (str:string) = 
    let tmp = str.Split([|'\t'|])
    {
    JGI43full           = tmp.[0]
    JGI55full           = tmp.[1]
    JGI43               = tmp.[2]
    JGI55               = tmp.[3]
    MapManGMM           = tmp.[4].Split([|','|])
    MapManDescription   = tmp.[5].Split([|','|])
    }

let readJGIMapping = 
    System.IO.File.ReadAllLines(@"c:\Users\mikha\Work-CSB\Data\jgi431_55_mapMan.txt")
    |> Array.tail
    |> Array.map (fun x -> 
        createJGIMapping x
        )

let mapFromHSREmicroArrayToMapManAnnotation =
    readJGIMapping
    |> Array.map (fun x -> x.JGI43,x.MapManGMM)
    |> Map.ofArray

//mapping from microarray cre number (jgi version 4.3) to MapMan Bin list (based on jgi version 5.5)
"Cre06.g250100" |> fun creNr -> if mapFromHSREmicroArrayToMapManAnnotation.ContainsKey creNr then mapFromHSREmicroArrayToMapManAnnotation.[creNr] else [|""|]


//type BinConverter() = 
//    inherit SchemaReader.Attribute.ConverterAttribute()
//    override this.convertToObj = 
//        SchemaReader.Converter.Single (fun (strs : string) -> 
//                                            (strs |> String.filter (fun c -> c<>''')|> String.split '.' ) |> box)

//type IdentifierConverter() = 
//    inherit SchemaReader.Attribute.ConverterAttribute()
//    override this.convertToObj = 
//        SchemaReader.Converter.Single (fun (strs : string) -> 
//            (strs 
//            |> String.split '|' 
//            |> Array.item 0 
//            |> String.replace "p" "" 
//            |> String.replace "cre" "Cre"
                                              
//            ) |> box)

type StringArrayConverter() = 
    inherit SchemaReader.Attribute.ConverterAttribute()
    override this.convertToObj = 
        SchemaReader.Converter.Collection (fun (strs : seq<string>) -> 
                                                (strs |> Seq.toArray) |> box)

// read MapMan to extract Bin

//type MapManRead = {
//    [<SchemaReader.Attribute.FieldAttribute("BINCODE")>]    [<BinConverter>]       BinString           : string []
//    [<SchemaReader.Attribute.FieldAttribute("IDENTIFIER")>] [<IdentifierConverter>]       ProteinIdentifier   : string
//    }

//let csvPathMM = sprintf @"%sdata\Creinhardtii_236 - ManMapList.txt" General.pathToData

//let readerMM    = new SchemaReader.Csv.CsvReader<MapManRead>(schemaMode=SchemaReader.Csv.SchemaMode.Fill)

//// Read Proteins database

type TranscriptItemReadShaped = {
    [<SchemaReader.Attribute.FieldAttribute("Description")>]                [<ChlamyProteome.NameConverter>]           ProteinGroup    : string []
    [<SchemaReader.Attribute.FieldAttribute(
        //[| "raw_means_0";"raw_means_1";"raw_means_2";"raw_means_3";
        //"raw_means_4";"raw_means_5";"raw_means_6";"raw_means_7"|])>]        [<ChlamyProteome.DoubleArrayConverter>]        Features        : float []
        [| "spline_A_0";"spline_A_1";"spline_A_2";"spline_A_3";
        "spline_A_4";"spline_A_5";"spline_A_6";"spline_A_7"|])>]        [<ChlamyProteome.DoubleArrayConverter>]        Features        : float []
    [<SchemaReader.Attribute.FieldAttribute([| "parent_shape";"extrema"|])>]    [<StringArrayConverter>]        ShapeClass        : string []
    }

let csvPathT = sprintf @"%sdata\HSRE_ChlamyTranscripts_with_derivativesUPDATE.txt" General.pathToData

let readerT    = new SchemaReader.Csv.CsvReader<TranscriptItemReadShaped>(schemaMode=SchemaReader.Csv.SchemaMode.Fill)

/// Variable, contains all raw data from protein DataBase file in csvPath
let dataProteinT = 
    readerT.ReadFile (csvPathT, General.separatorTab, General.hasHeader) 
    |> Seq.distinctBy (fun x -> x.ProteinGroup.[0])
    |> Seq.toArray

dataProteinT.Length

let hermitShapeMap =
    let temp =
        dataProteinT
        |> Array.map (fun x -> 
            let shape = 
                if x.ShapeClass.[0]="In0" then [|(-1.,-1.)|]
                elif x.ShapeClass.[0]="De0" then [|(1.,-1.)|]
                else
                    x.ShapeClass.[1] 
                    |> String.replace "\"" ""
                    |> String.replace "[" "" 
                    |> String.replace "]" "" 
                    |> String.split ';'
                    //|> Array.filter (fun point -> point |> String.contains ",")
                    |> Array.map (fun x -> 
                        x
                        |> String.replace "(" "" 
                        |> String.replace ")" "" 
                        |> String.replace "," "." 
                        |> String.toWords
                        |> Seq.toArray
                        |> (fun x -> x.[0] |> String.tryParseFloatDefault 0. , x.[1] |> String.tryParseFloatDefault 0. |> round 0 )
                        )
            x.ProteinGroup.[0], shape)
    let patterns = temp |> Array.distinctBy (snd) |> Array.mapi (fun id (_,s) -> (id,s))
    temp
    |> Array.map (fun (p,s) -> (p, patterns |> Array.find (fun (id,sid) -> sid=s) |> fst))
    |> Map.ofArray

///// Data to find mapman bin for proteins
//let dataMM = 
//    readerMM.ReadFile (csvPathMM, General.separatorTab, General.hasHeader) 
//    |> Seq.filter (fun x -> x.ProteinIdentifier.Length > 10 && x.ProteinIdentifier.Contains ".")
//    |> Seq.distinctBy (fun x -> x.ProteinIdentifier)
//    |> Seq.map (fun x -> {x with ProteinIdentifier = x.ProteinIdentifier |> String.subString 1 13  })
//    |> Seq.toArray

let dataWithBin =
    dataProteinT
    |> Array.choose (fun i -> 
        try Some 
                    
                (
                //let bin = (dataMM |> Array.find (fun ii -> ii.ProteinIdentifier = i.ProteinGroup.[0]) ).BinString
                let bin = 
                    (readJGIMapping |> Array.find (fun ii -> ii.JGI43 = i.ProteinGroup.[0]) ).MapManGMM.[0]
                    |> String.split ':'
                    |> Array.item 1
                    |> String.split '.'

                {
                ID = 0;
                ProteinL = i.ProteinGroup;
                OriginalBin = bin;
                BinL = bin;
                dataL = General.zScoreTransform i.Features 
                } )

                with

            |ex -> None)

dataProteinT.Length
dataWithBin.Length // 9330 vs 5580

let dataGroupedHermit =
    dataWithBin 
    |> Array.groupBy (fun x -> Map.find x.ProteinL.[0] hermitShapeMap)

dataGroupedHermit
|> Array.map (fun i -> i |> snd |> Array.length)

//
////---------------------------------------------Insert Hermit group in MapMan raw--------------------------------------------------------------------
//

let rec breakGroupHermit (hermitShapeMap: Map<string,int>) (items: Types.Item array) depth =
    
    //printfn "break for depth %i" depth
    //printfn "%A" (items |> Array.map (fun i -> i.BinL))

    if (items |> Array.forall (fun ii -> ii.OriginalBin.Length>depth+1 && ii.OriginalBin.[depth+1]=items.[0].OriginalBin.[depth+1])) 
    then // there is a tunnel -> delete the tunneling subbin
        
        let newItems = 
            items
            |> Array.map (fun x -> {x with 
                                        OriginalBin=(x.OriginalBin |> Array.removeIndex (depth))
                                        BinL=(x.BinL |> Array.removeIndex (depth))})
        breakGroupHermit hermitShapeMap newItems depth

    elif (depth>0 && items.[0].BinL.[depth] |> String.contains "h" )
    then // we are inside Hermit group -> break into singletons

        items 
        |> Array.map (fun i -> 
            if i.BinL.Length<=(depth+1) then 
                {i with 
                    OriginalBin = Array.append i.OriginalBin [|sprintf "%s%i" "p" i.ID|]; 
                    BinL = Array.append i.BinL [|sprintf "%s%i" "p" i.ID|]}
            else i)
        |> Array.groupBy (fun i -> i.BinL.[depth+1])
        |> Map.ofArray

    else // no tunnels in the structure -> make a break

        items 
        |> Array.map (fun i -> 
            if i.BinL.Length<=(depth+1) then 
                let hermitClass = hermitShapeMap |> Map.find i.ProteinL.[0]
                {i with 
                    //OriginalBin = Array.append i.OriginalBin [|sprintf "%s%i" "h" hermitClass|]; 
                    BinL = Array.append i.BinL [|sprintf "%s%i" "h" hermitClass|]}
            else i)
        |> Array.groupBy (fun i -> i.BinL.[depth+1])
        |> Map.ofArray

let createTree gainFn (weight: seq<float> option) (mode: int) (rootGroup: Types.Item array) = 
        
    let nRoot = rootGroup.Length

    let matrix = 
        rootGroup
        |> General.distMatrixWeightedOf General.distanceMatrixWeighted weight

    // calculation for one node    
    let rec loop (nodeMembers: Types.Item array) depth dPredSum =
        
        /// sum of max dist within a node 
        let dCurrSum = General.dSumFn (General.groupIDFn nodeMembers) (General.groupIDFn nodeMembers) matrix

        /// to calc step from parent to current node
        let stepGain = (gainFn dCurrSum dPredSum nodeMembers.Length nRoot)

        let children = 
            match mode with
            |0 -> // mapman with Hermit layer
                if nodeMembers.Length=1 then
                    Map.empty
                else 
                    (breakGroupHermit hermitShapeMap nodeMembers depth)
                    |> Map.map (fun key nodes -> 
                                    let dPredSum' = General.dSumFn (General.groupIDFn nodes) (General.groupIDFn nodeMembers) matrix
                                    (loop nodes (depth+1) dPredSum'))


            |1 -> // combi without Hermit layer
                if (nodeMembers.Length=1) 
                        || (nodeMembers.Length=0) 
                        || (String.contains "mix" (String.Concat nodeMembers.[0].BinL)) 
                        || (String.contains "c" (String.Concat nodeMembers.[0].BinL)) 
                        || (String.contains "|" (String.Concat nodeMembers.[0].BinL)) then 
                    Map.empty
                else 
                    (General.breakGroup nodeMembers depth)
                    |> General.partGroup depth
                    |> Seq.fold (fun (singles,best) i -> 
                        let newNodes = 
                            if singles=Map.empty then
                                i
                                |> Map.fold (fun state key nodes ->
                                    let dPredSum' = General.dSumFn (General.groupIDFn nodes) (General.groupIDFn nodeMembers) matrix
                                    state 
                                    |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)
                            else
                                i
                                |> Map.fold (fun state key nodes ->
                                    match (singles.TryFind key) with
                                    | None ->
                                        let dPredSum' = General.dSumFn (General.groupIDFn nodes) (General.groupIDFn nodeMembers) matrix
                                        state 
                                        |> Map.add key (loop nodes (depth+1) dPredSum')
                                    | Some x ->
                                        state 
                                        |> Map.add key x                                                                   
                                ) (Map.empty)
                                                                    
                        let best' =
                            if (General.confGainFn newNodes) > (General.confGainFn best) then  // compare configuration gains to get the best
                                newNodes  
                            else 
                                best
                        if (singles = Map.empty) then
                            (newNodes, best')
                        else
                            (singles, best')
                                            
                    ) (Map.empty,Map.empty) // here as state should be this singles (first) saved and optimal conf
                    |> snd

            |2 -> // combi with Hermit layer
                if (nodeMembers.Length=1) 
                        || (nodeMembers.Length=0) 
                        || (String.contains "mix" (String.Concat nodeMembers.[0].BinL)) 
                        || (String.contains "c" (String.Concat nodeMembers.[0].BinL)) 
                        || (String.contains "|" (String.Concat nodeMembers.[0].BinL)) then 
                    Map.empty
                else 
                    (breakGroupHermit hermitShapeMap nodeMembers depth)
                    |> General.partGroup depth
                    |> Seq.fold (fun (singles,best) i -> 
                        let newNodes = 
                            if singles=Map.empty then
                                i
                                |> Map.fold (fun state key nodes ->
                                    let dPredSum' = General.dSumFn (General.groupIDFn nodes) (General.groupIDFn nodeMembers) matrix
                                    state 
                                    |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)
                            else
                                i
                                |> Map.fold (fun state key nodes ->
                                    match (singles.TryFind key) with
                                    | None ->
                                        let dPredSum' = General.dSumFn (General.groupIDFn nodes) (General.groupIDFn nodeMembers) matrix
                                        state 
                                        |> Map.add key (loop nodes (depth+1) dPredSum')
                                    | Some x ->
                                        state 
                                        |> Map.add key x                                                                   
                                ) (Map.empty)
                                                                    
                        let best' =
                            if (General.confGainFn newNodes) > (General.confGainFn best) then  // compare configuration gains to get the best
                                newNodes  
                            else 
                                best
                        if (singles = Map.empty) then
                            (newNodes, best')
                        else
                            (singles, best')
                                            
                    ) (Map.empty,Map.empty) // here as state should be this singles (first) saved and optimal conf
                    |> snd

            |_ -> // with Hier clust as a start point for gain walking and with Hermit layer
                if (nodeMembers.Length=1) 
                        || (nodeMembers.Length=0) 
                        || (String.contains "|" (String.Concat nodeMembers.[0].BinL)) then 
                    Map.empty
                else 
                    (breakGroupHermit hermitShapeMap nodeMembers depth)
                    |> (fun x -> 
                                                
                        let singles = 
                            x
                            |> Map.fold (fun state key nodes ->
                                    let dPredSum' = General.dSumFn (General.groupIDFn nodes) (General.groupIDFn nodeMembers) matrix
                                    state 
                                    |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)

                        (Walk.walkingHC_Fn Clustering.clusterHierGroups) depth matrix singles gainFn x 
                        |> Map.fold (fun state key nodes ->
                            match (singles.TryFind key) with
                            | None ->
                                let dPredSum' = General.dSumFn (General.groupIDFn nodes) (General.groupIDFn nodeMembers) matrix
                                state 
                                |> Map.add key (loop nodes (depth+1) dPredSum')
                            | Some x ->
                                state 
                                |> Map.add key x                                                                   
                        ) (Map.empty) )

        let confGain = General.confGainFn children
        {
        Member = nodeMembers;
        Children = 
            match mode with
            |0 -> 
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

let createTreeVirtualNodes gainFn (weight: seq<float> option) (mode: int) (rootGroup: Types.Item array) = 
        
    let nRoot = rootGroup.Length

    let matrix = 
        rootGroup
        |> General.distMatrixWeightedOf General.distanceMatrixWeighted weight

    // calculation for one node    
    let rec loop (nodeMembers: Types.Item array) depth dPredSum =
        
        /// sum of max dist within a node 
        let dCurrSum = General.dSumFn (General.groupIDFn nodeMembers) (General.groupIDFn nodeMembers) matrix

        /// to calc step from parent to current node
        let stepGain = (gainFn dCurrSum dPredSum nodeMembers.Length nRoot)

        let children = 
            match mode with
            |0 -> // mapman with Hermit layer
                if nodeMembers.Length=1 then
                    Map.empty
                else 
                    (breakGroupHermit hermitShapeMap nodeMembers depth)
                    |> Map.map (fun key nodes -> 
                                    let dPredSum' = General.dSumFn (General.groupIDFn nodes) (General.groupIDFn nodeMembers) matrix
                                    (loop nodes (depth+1) dPredSum'))


            |1 -> // combi without Hermit layer
                if (nodeMembers.Length=1) 
                        || (nodeMembers.Length=0) 
                        || (String.contains "mix" (String.Concat nodeMembers.[0].BinL)) 
                        || (String.contains "c" (String.Concat nodeMembers.[0].BinL)) 
                        || (String.contains "|" (String.Concat nodeMembers.[0].BinL)) then 
                    Map.empty
                else 
                    (General.breakGroup nodeMembers depth)
                    |> General.partGroup depth
                    |> Seq.fold (fun (singles,best) i -> 
                        let newNodes = 
                            if singles=Map.empty then
                                i
                                |> Map.fold (fun state key nodes ->
                                    let dPredSum' = General.dSumFn (General.groupIDFn nodes) (General.groupIDFn nodeMembers) matrix
                                    state 
                                    |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)
                            else
                                i
                                |> Map.fold (fun state key nodes ->
                                    match (singles.TryFind key) with
                                    | None ->
                                        let dPredSum' = General.dSumFn (General.groupIDFn nodes) (General.groupIDFn nodeMembers) matrix
                                        state 
                                        |> Map.add key (loop nodes (depth+1) dPredSum')
                                    | Some x ->
                                        state 
                                        |> Map.add key x                                                                   
                                ) (Map.empty)
                                                                    
                        let best' =
                            if (General.confGainFn newNodes) > (General.confGainFn best) then  // compare configuration gains to get the best
                                newNodes  
                            else 
                                best
                        if (singles = Map.empty) then
                            (newNodes, best')
                        else
                            (singles, best')
                                            
                    ) (Map.empty,Map.empty) // here as state should be this singles (first) saved and optimal conf
                    |> snd

            |2 -> // combi with Hermit layer
                if (nodeMembers.Length=1) 
                        || (nodeMembers.Length=0) 
                        || (String.contains "mix" (String.Concat nodeMembers.[0].BinL)) 
                        || (String.contains "c" (String.Concat nodeMembers.[0].BinL)) 
                        || (String.contains "|" (String.Concat nodeMembers.[0].BinL)) then 
                    Map.empty
                else 
                    (breakGroupHermit hermitShapeMap nodeMembers depth)
                    |> General.partGroup depth
                    |> Seq.fold (fun (singles,best) i -> 
                        let newNodes = 
                            if singles=Map.empty then
                                i
                                |> Map.fold (fun state key nodes ->
                                    let dPredSum' = General.dSumFn (General.groupIDFn nodes) (General.groupIDFn nodeMembers) matrix
                                    state 
                                    |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)
                            else
                                i
                                |> Map.fold (fun state key nodes ->
                                    match (singles.TryFind key) with
                                    | None ->
                                        let dPredSum' = General.dSumFn (General.groupIDFn nodes) (General.groupIDFn nodeMembers) matrix
                                        state 
                                        |> Map.add key (loop nodes (depth+1) dPredSum')
                                    | Some x ->
                                        state 
                                        |> Map.add key x                                                                   
                                ) (Map.empty)
                                                                    
                        let best' =
                            if (General.confGainFn newNodes) > (General.confGainFn best) then  // compare configuration gains to get the best
                                newNodes  
                            else 
                                best
                        if (singles = Map.empty) then
                            (newNodes, best')
                        else
                            (singles, best')
                                            
                    ) (Map.empty,Map.empty) // here as state should be this singles (first) saved and optimal conf
                    |> snd

            |_ -> // with Hier clust as a start point for gain walking and with Hermit layer
                if (nodeMembers.Length=1) 
                        || (nodeMembers.Length=0) 
                        || (String.contains "|" (String.Concat nodeMembers.[0].BinL)) then 
                    Map.empty
                else 
                    (breakGroupHermit hermitShapeMap nodeMembers depth)
                    |> (fun x -> 
                                                
                        let singles = 
                            x
                            |> Map.fold (fun state key nodes ->
                                    let dPredSum' = General.dSumFn (General.groupIDFn nodes) (General.groupIDFn nodeMembers) matrix
                                    state 
                                    |> Map.add key (loop nodes (depth+1) dPredSum')) (Map.empty)

                        (Walk.walkingHC_Fn Clustering.clusterHierGroups) depth matrix singles gainFn x 
                        |> Map.fold (fun state key nodes ->
                            match (singles.TryFind key) with
                            | None ->
                                let dPredSum' = General.dSumFn (General.groupIDFn nodes) (General.groupIDFn nodeMembers) matrix
                                state 
                                |> Map.add key (loop nodes (depth+1) dPredSum')
                            | Some x ->
                                state 
                                |> Map.add key x                                                                   
                        ) (Map.empty) )

        let clearHfromSubbin (node: Node<string,Item>) =
            node.Member
            |> Array.map (fun x -> 
                if x.OriginalBin.Length<=(depth+1) then 
                    sprintf "p%i" x.ID
                else
                    x.OriginalBin.[depth+1])
            |> Array.distinct
            |> fun ar ->
                if ar.Length=1 then 
                    sprintf "%s" ar.[0]
                else
                    ar
                    |> Array.fold (fun acc x -> sprintf "%s|%s" acc x) "mix"

        let changeLastSubbin newBin (node: Node<string,Item>) =
            let lastN = node.Member.[0].BinL.Length - 1
            let x = Array.copy node.Member.[0].BinL
            x.[lastN] <- newBin
            x

        let removeBinBeforeLast item =
            let beforeLastN = item.BinL.Length - 2
            let x = Array.copy item.BinL
            x |> Array.removeIndex beforeLastN

        let children_noHermit =
            children
            |> Map.toArray
            |> Array.map (fun (bin,node) -> 
                if (bin |> String.contains "h") then
                    if node.Children=Map.empty then
                        //only delete h from last bin name in bin and members
                        let newBin = clearHfromSubbin node
                        let newBinL = changeLastSubbin newBin node
                        [|(newBin,{node with Member=(node.Member |> Array.map (fun item -> {item with BinL=newBinL}))})|]
                    else
                        // move child's children one level higher with renaming members' bins and recalculating step gains
                        (node.Children 
                        |> Map.toArray 
                        |> Array.map (fun (subbin,subnode) -> 
                            let newBin = removeBinBeforeLast subnode.Member.[0]
                            let newStepGain = General.getStepGainFn gainFn (General.groupIDFn subnode.Member) (General.groupIDFn nodeMembers) nodeMembers.Length matrix
                            (subbin, {subnode with 
                                        Member=(subnode.Member |> Array.map (fun subitem -> {subitem with BinL=newBin})); 
                                        StepGain=newStepGain;
                                        GroupGain=newStepGain})))
                else
                    // no changes
                    [|(bin,node)|] )
            |> Array.concat
            |> Map.ofArray

        let confGain = General.confGainFn children_noHermit
        {
        Member = nodeMembers;
        Children = 
            match mode with
            |0 -> 
                children_noHermit
            |_ -> 
                if  (confGain > stepGain) then 
                    children_noHermit;
                else 
                    Map.empty
        StepGain = stepGain; 
        ConfGain = (confGain, children_noHermit |> Map.toList |> List.map fst);
        GroupGain = max stepGain confGain;
        }
    
    loop rootGroup 0 0.

//
////----------------------------------------Test-----------------------------------------------------------------------------------------------------------
//

let path_of_interest =
    dataWithBin 
    |> Array.filter (fun x -> x.BinL.[0]="13")
    |> Array.mapi (fun id x -> {x with ID=id})

path_of_interest.Length

path_of_interest |> Array.find (fun i -> i.ProteinL |> Array.contains "Cre05.g237450")

let mapman29 = readMM path_of_interest.Length path_of_interest
let mapman29_R = readMM_raw path_of_interest.Length path_of_interest


mapman29 |> sendToGephiFromTreeParam
mapman29_R |> sendToGephiFromTreeParam

mapman29 |> Tree.filterChildrenList |> List.sortDescending

let biggestNodeSingle = mapman29 |> Tree.findNodeMembers ["29";"4"] |> Array.filter (fun i -> i.BinL.Length=2)
let biggestNodeSingleHermit =
    biggestNodeSingle |> Array.groupBy (fun x -> Map.find x.ProteinL.[0] hermitShapeMap)

biggestNodeSingle.Length
biggestNodeSingleHermit |> Array.map (snd >> Array.length)


let sstHermit_29 = createTree (getStepGainNodeSetnR path_of_interest.Length)  (None) 10 path_of_interest

let sstCombi_29 = createTree (getStepGainNodeSetnR path_of_interest.Length)  (None) 1 path_of_interest

sstHermit_29.GroupGain
sstHermit_29 |> sendToGephiFromTreeParam

sstHermit_29 |> Tree.filterChildrenList |> List.sortDescending

drawKinetikTitle sstHermit_29 ["29";"2";"1";"1";"1";"1"] [|0 .. 7|] |> Chart.Show

let node = sstHermit_29 |> Tree.findNode ["29";"2";"1";"1";"1";"1"] |> Tree.filterLeaves
drawKinetikRange [|0. .. 7.|] "Ribosomal Proteins vs PSRP1 in [29;2;1;1;1;1]" node |> Chart.Show
drawLeaves "Ribosomal Proteins vs PSRP1 in [29;2;1;1;1;1]" [|0. .. 7.|] (sstHermit_29 |> Tree.findNode ["29";"2";"1";"1";"1";"1"])

sstCombi_29.GroupGain
sstCombi_29 |> sendToGephiFromTreeParam

Tree.treeComparison sstHermit_29 sstCombi_29

let mmH_poi = (createTree (getStepGainNodeSetnR path_of_interest.Length) (None) 0 path_of_interest)

let combi_poi = (createTree (getStepGainNodeSetnR path_of_interest.Length) (None) 1 path_of_interest)
let combiH_poi = (createTree (getStepGainNodeSetnR path_of_interest.Length) (None) 2 path_of_interest)

let combiH_poi_virtual = (createTreeVirtualNodes (getStepGainNodeSetnR path_of_interest.Length) (None) 2 path_of_interest)

let sstHvirtual_poi_ = (createTreeVirtualNodes (getStepGainNodeSetnR path_of_interest.Length) (None) 10 path_of_interest)
let sstHermit_poi = createTree (getStepGainNodeSetnR path_of_interest.Length)  (None) 10 path_of_interest
sstHermit_poi.GroupGain
sstHvirtual_poi_.GroupGain
sstHermit_poi |> sendToGephiFromTreeParam
sstHvirtual_poi_ |> sendToGephiFromTreeParam
sstHvirtual_poi_ |> Tree.findNode ["13";"2";"4";"1"]
combi_poi |> Tree.findNode ["13";"2";"4";"1"]
mmH_poi |> Tree.findNode ["13";"2";"4";"1"]

sstHvirtual_poi_ |> Tree.filterLeaves |> Array.concat |> Array.length

drawLeaves "path 13, trouble subbin [13;2;4;1]" [|0. .. 7.|] (mmH_poi |> Tree.findNode ["13";"2";"4";"1"])

hermitShapeMap |> Map.find path_of_interest.[8].ProteinL.[0]  // Cre13.g596350  1   De0	[]
hermitShapeMap |> Map.find path_of_interest.[39].ProteinL.[0] // Cre18.g746550  13  In2	[(1, 1,84005872078159); (-1, 5,01705742658671)]
hermitShapeMap |> Map.find path_of_interest.[23].ProteinL.[0] // Cre13.g596250  3   In2	[(1, 1,32341395811728); (-1, 4,99999999999915)]
hermitShapeMap |> Map.find path_of_interest.[21].ProteinL.[0] // Cre13.g596650  3   In2	[(1, 1); (-1, 2,41711557416432)]
hermitShapeMap |> Map.find path_of_interest.[30].ProteinL.[0] // Cre18.g746750  3   In2	[(1, 1); (-1, 2,99999998376502)]

let check = dataProteinT |> Array.find (fun i -> i.ProteinGroup.[0]="Cre18.g746550")

"1.84005872078159" |> String.tryParseFloatDefault 0. |> round 0

check.ShapeClass.[1]
|> String.replace "[" "" 
|> String.replace "]" "" 
|> String.split ';'
|> Array.map (fun x -> 
    x
    |> String.replace "\"" ""
    |> String.replace "(" "" 
    |> String.replace ")" "" 
    |> String.replace "," "." 
    |> String.toWords
    |> Seq.toArray
    |> (fun x -> x.[0] |> String.tryParseFloatDefault 0. , x.[1] |> String.tryParseFloatDefault 0. |> round 0 )
    )

[[|(1.0, 2.0); (-1.0, 5.0)|];[|(1.0, 1.0); (-1.0, 5.0)|];[|(1.0, 1.0); (-1.0, 2.0)|];[|(1.0, 1.0); (-1.0, 3.0)|]] |> List.distinct

combi_poi |> sendToGephiFromTreeParam
combiH_poi |> sendToGephiFromTreeParam

combiH_poi_virtual |> sendToGephiFromTreeParam


//
////--------------------------------------------Analysis-----------------------------------------------------------------
//

let pathsSorted =
    dataWithBin
    |> Array.groupBy (fun x -> x.BinL.[0])
    |> Array.filter (fun (bin,l) -> l.Length > 2 && bin<>"35")
    |> Array.map (fun (bin,l) -> 
        let data = l |> Array.mapi (fun id x -> {x with ID=id})
        let tree = readMM data.Length data
        let treeMMHermit = createTree (getStepGainNodeSetnR data.Length) (None) 0 data
        let childrenN = Tree.filterChildrenList tree |> List.max
        let childrenNHermit = Tree.filterChildrenList treeMMHermit |> List.max

        let (combiHermitGG, cHTime, cHdxc) =
            if (childrenNHermit)>10 then 
                0., 0., 0.
            else
                let stopwatch = new System.Diagnostics.Stopwatch()
                stopwatch.Start()
                let combiTree = (createTreeVirtualNodes (getStepGainNodeSetnR data.Length) (None) 2 data) // createTree (getStepGainNodeSetnR data.Length) (None) 1 data
                let combiTime = (stopwatch.Elapsed.TotalSeconds)
                stopwatch.Stop()

                let matrixOfN = 
                    data
                    |> distMatrixWeightedOf distanceMatrixWeighted None
            
                let dissMax = 
                    matrixOfN
                    |> Array2D.array2D_to_seq
                    |> Seq.max

                let normMatrix = 
                    matrixOfN
                    |> Array2D.map (fun i -> i/dissMax)

                let combidxc = 
                    combiTree
                    |> Tree.filterLeaves
                    |> Analysis.pointDxC normMatrix
                    |> fun (x,y) -> weightedEuclidean None [x;y] [0.;0.]
                (combiTree.GroupGain, combiTime, combidxc)
                
                //(createTree (getStepGainNodeSetnR data.Length) (None) 2 data).GroupGain

        printfn "%s\t%i\t%i\t%f\t%i\t%f\t%f\t%f\t%f" bin data.Length childrenN tree.GroupGain childrenNHermit treeMMHermit.GroupGain combiHermitGG cHTime cHdxc
        (bin, data.Length, childrenN, tree.GroupGain, childrenNHermit, treeMMHermit.GroupGain ) )
    


let pathFile = sprintf @"%sHermitLayer\" @"c:\Users\mikha\Work-CSB\Projects\SSN\results\"
let neader = "Path\tnRoot\tchildrenMax\tcombi_GG\tcombi_Time\tcombi_DxC\twalk_GG\twalk_Time\twalk_DxC"

File.AppendAllLines((sprintf "%s%s.txt" pathFile "ChamyTranscrData_combi_vs_walk_d1_best1_del06"), [neader])

let linesWAcombi = 
    dataWithBin
    |> Array.groupBy (fun x -> x.BinL.[0])
    |> Array.filter (fun (bin,l) -> l.Length > 2 && bin<>"35")
    |> Array.map (fun (bin,l) -> 
        let data = 
            l
            |> Array.mapi (fun id x -> {x with ID=id})
        let tree = readMM data.Length data
        let childrenN = Tree.filterChildrenList tree |> List.max
        let treeMMHermit = createTree (getStepGainNodeSetnR data.Length) (None) 0 data
        let childrenNHermit = Tree.filterChildrenList treeMMHermit |> List.max
        let stopwatch = new System.Diagnostics.Stopwatch()
        
        let (combiGG,combiTime,combiDxC) = 
            if childrenNHermit>10 then (0.,0.,0.)
            else
                printfn "start combi for path %s" bin
                stopwatch.Start()
                let combiTree = (createTree (getStepGainNodeSetnR data.Length) (None) 2 data) // createTree (getStepGainNodeSetnR data.Length) (None) 1 data
                let combiTime = (stopwatch.Elapsed.TotalSeconds)
                stopwatch.Stop()

                let matrixOfN = 
                    data
                    |> distMatrixWeightedOf distanceMatrixWeighted None
            
                let dissMax = 
                    matrixOfN
                    |> Array2D.array2D_to_seq
                    |> Seq.max

                let normMatrix = 
                    matrixOfN
                    |> Array2D.map (fun i -> i/dissMax)

                let combidxc = 
                    combiTree
                    |> Tree.filterLeaves
                    |> Analysis.pointDxC normMatrix
                    |> fun (x,y) -> weightedEuclidean None [x;y] [0.;0.]
                printfn "end combi for path %s" bin
                (combiTree.GroupGain, combiTime, combidxc)

        printfn "start walkHermit for path %s" bin
        stopwatch.Start()
        let walk = createTreeVirtualNodes (getStepGainNodeSetnR data.Length) (None) 10 data // createTree (getStepGainNodeSetnR data.Length) (None) 10 data // 
        let timeWalk = (stopwatch.Elapsed.TotalSeconds)
        stopwatch.Stop()

        let matrixOfN = 
            data
            |> distMatrixWeightedOf distanceMatrixWeighted None
            
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
        
        printfn "path; nRoot; childrenMax; combi GG; combi Time; walk GG; walk Time; walk DxC"
        printfn "%s %i %i %f %f %f %f %f %f" 
            bin data.Length childrenNHermit combiGG combiTime combiDxC walk.GroupGain timeWalk dxcW

        let treeLinesW =
            walk
            |> Tree.filterLeaves
            |> Array.concat
            |> Array.map (fun x ->  (sprintf "%i" x.ID) + "\t" + (sprintf "%A" x.ProteinL) + "\t" + (String.Join(";", x.OriginalBin))  + "\t" + (String.Join(";", x.BinL)))
        File.AppendAllLines(sprintf @"c:\Users\mikha\Work-CSB\Projects\SSN\results\HermitLayer\ChlamyTranscr_HermitWalk_1d006del1best_%s.txt" bin, treeLinesW)


        let line = 
            sprintf "%s\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f" 
                bin data.Length childrenNHermit combiGG combiTime combiDxC walk.GroupGain timeWalk dxcW
        
        File.AppendAllLines((sprintf "%s%s.txt" pathFile "ChamyTranscrData_combi_vs_walk_d1_best1_del06"), [line])

        line
        )

