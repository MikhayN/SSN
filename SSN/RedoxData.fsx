#r @".\bin\Debug\netstandard.dll"

#r @".\bin\Debug\BioFSharp.dll"
#r @".\bin\Debug\BioFSharp.IO.dll"
#r @".\bin\Debug\FSharpAux.dll"
#r @".\bin\Debug\FSharpAux.IO.dll"
#r @".\bin\Debug\Newtonsoft.Json.dll"
#r @".\bin\Debug\FSharp.Stats.dll"
#r @".\bin\Debug\FSharp.Plotly.dll"
#r @".\bin\Debug\FSharpGephiStreamer.dll"

#r @"..\lib\MathNet.Numerics.dll"

#time

open System
open System.IO
open FSharp.Stats
open FSharp.Plotly
open BioFSharp.IO

open FSharpAux
open FSharpAux.IO


open BioFSharp
open BioFSharp
open BioFSharp.Elements

let findMatches chars str =
    /// Returns whether or not the string matches the beginning of the character array
    let rec isStartMatch chars (str: string) =
        match chars with
        | char :: rest when str.Length > 0 ->
            char = str.[0] && (isStartMatch rest str.[1..(str.Length - 1)])
        | _ -> str.Length = 0
    /// The actual function here
    let rec findMatches matchedIndices i chars str =
        match chars with
        | _ :: rest ->
            if isStartMatch chars str
            then findMatches (i :: matchedIndices) (i + 1) rest str
            else findMatches matchedIndices (i + 1) rest str
        | [] -> matchedIndices

    findMatches [] 0 chars str

let adjustCysPosition protSeq cysSeq cysPs a =
    if protSeq |> String.contains cysSeq then
        let c = cysPs |> Array.findIndex (fun i -> i=a)
        let pInSub = cysSeq |> String.toCharArray |> Array.indexed |> Array.filter (fun (i,c) -> c='C') |> Array.item c |> fst
        let startSubseq = findMatches (protSeq |> String.toCharArray |> Array.toList) cysSeq
        startSubseq |> List.map (fun x -> x + pInSub) |> List.minBy (fun t -> abs (t-a))
    else a

//Read file of RedoxDB.A

type CysItem =
    {
    CysID: int;
    Organism: string;
    ProtID: string;
    ProteinSeq: string;
    CysPosition: int;
    CysSeq: string;
    Label: bool;
    }

let fillCysItem cysID org protID protSeq cysP cysSeq label = {CysID=cysID; Organism=org; ProtID=protID; ProteinSeq=protSeq; CysPosition=cysP; CysSeq=cysSeq; Label=label}

let enum id item = {item with CysID=id}

type FastaItem = {
    Header    : string;
    Sequence  : string;       
    }

/// Creates with header line and sequence.
let createFastaItem header sequence =
    { Header = header; Sequence = sequence }

// Conditon of grouping lines
let same_group l =             
    not (String.length l = 0 || l.[0] <> '>')
    
/// Reads FastaItem from file. Converter determines type of sequence by converting seq<char> -> type
let fromFileEnumerator (converter:seq<char>-> string) (fileEnumerator) =
    // Matches grouped lines and concatenates them
    let record d (converter:seq<char>-> string) = 
        match d with
        | [] -> raise (System.Exception "Incorrect FASTA format")
        | (h:string) :: t when h.StartsWith ">" ->  let header = h .Remove(0,1)
                                                    let sequence = (Seq.concat t) |> converter
                                                    createFastaItem header sequence
                                                        
        | h :: _ -> raise (System.Exception "Incorrect FASTA format")        
    
    // main
    fileEnumerator
    |> Seq.filter (fun (l:string) -> not (l.StartsWith ";" || l.StartsWith "#"))
    |> Seq.groupWhen same_group 
    |> Seq.map (fun l -> record (List.ofSeq l) converter)

/// Reads FastaItem from file. Converter determines type of sequence by converting seq<char> -> type
let fromFile converter (filePath) =
    FileIO.readFile filePath
    |> fromFileEnumerator converter


////////////// read RedoxData

let linesA = FileIO.readFile @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\redoxdbA.txt"

let itemsA =
    linesA 
    |> Seq.fold (fun (item,list) i -> 
        if i="" then ([],item::list)
        else (i::item,list)) ([],[])
    |> snd

let readA (line: string list) =
    let org = line |> List.find (fun s -> s |> String.startsWith "ORGANISM") |> String.split ':' |> Array.item 1
    let protID = line.[1] |>  String.replace ">" "" 
    let protSeq = if (line.[0] |> String.startsWith protID) then (line.[0].[protID.Length ..]) else line.[0]
    let cysP = line |> List.filter (fun s -> s |> String.startsWith "CYSTEINE") |> List.map (String.split ' ' >> Array.item 1 >> String.replace "," "" >> String.tryParseIntDefault 0)
    let label = true
    [for a in cysP -> 
        let cP = 
            if protSeq.[a-1]='C' then a-1
            elif protSeq.[a]='C' then a
            elif protSeq.[a+1]='C' then a+1
            elif protSeq.[a+2]='C' then a+2
            elif protSeq.[a-2]='C' then a-2
            elif protSeq.[a+3]='C' then a+3
            elif protSeq.[a-3]='C' then a-3
            else a
        printfn "seq_Length=%i, item at %i is %c" protSeq.Length cP protSeq.[cP]
        fillCysItem 0 org protID protSeq (cP) "" label]

let redoxPart =
    itemsA
    |> List.filter (fun x -> x<>[])
    |> List.map readA
    |> List.concat
    |> List.toArray
    |> Array.filter (fun i -> i.ProteinSeq.[i.CysPosition]='C')


let findNegative (protSeq: string) cysPP = ///// check the function, if the CysPosition is coreecly indexed!
    let cysPs = protSeq |> String.toCharArray |> Array.indexed |> Array.filter (fun (i,res) -> res='C') |> Array.map (fun (x,_) -> x+1)
    cysPP |> Set.ofArray |> Set.difference (cysPs |> Set.ofArray) |> Set.toArray

let createNegative (dbPart: CysItem []) (seqItem: CysItem)  =
    let cysPP = dbPart |> Array.filter (fun i -> i.ProtID = seqItem.ProtID ) |> Array.map (fun i -> i.CysPosition)
    let cysPN = findNegative seqItem.ProteinSeq cysPP
    cysPN
    |> Array.map (fun i -> fillCysItem 0 seqItem.Organism seqItem.ProtID seqItem.ProteinSeq i "" false)

let negativePool = 
    redoxPart
    |> Array.distinctBy (fun i -> i.ProtID)
    |> Array.map (createNegative redoxPart)
    |> Array.concat

let negativeRedoxPart =
    let r = new Random()
    Array.sampleWithOutReplacement r (negativePool) (redoxPart.Length)

let redoxData =
    [redoxPart;negativeRedoxPart]
    |> Array.concat
    |> Array.mapi enum
redoxPart.Length

/// read Lemair data

let fileDirLem = @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\Lemaire_2017.fasta"

//reads from file to an array of FastaItems.
let sequencesLem = 
    fileDirLem
    |> fromFile (Seq.toArray >> String.fromCharArray)
    |> Seq.toArray
    |> Array.map (fun i -> 
            let h = i.Header |> String.split '|' |> Array.item 1 //|> String.tryParseIntDefault 0
            {i with Header=h})

let linesC1 = FileIO.readFile @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\Lemaire_2017.txt"

let lines =
    linesC1
    |> Seq.toList  
    |> List.map (fun x -> (x |> String.split '\t'))
let llines = 
    lines.[5 ..] 
    |> List.map (fun i -> [|i.[14];i.[6];i.[5]|]) 
    |> List.filter (fun i -> i.[1]<>"")
let lllines =
    llines
    |> List.mapFold (fun lastProtID i -> 
        if i.[0]="" then
            [|lastProtID;i.[1];i.[2]|], lastProtID
        else
            i, i.[0]
        ) ""
    |> fst

let itemsC =
    lllines
    |> List.map (fun x -> 
        let org = "Chlamy"
        let protID = x.[0]
        let protSeq = (sequencesLem |> Array.find (fun i -> i.Header=protID)).Sequence
        let cysSeq = x.[2]
        let cysPs = 
            x.[1] 
            |> String.replace " ; " ";" 
            |> String.replace "C" "" 
            |> String.replace "\"" ""
            |> String.split ';'  
            |> Array.map (String.tryParseIntDefault 0)
        let label = true
        [for a in cysPs -> 
            let cP = 
                if a>=protSeq.Length-1  then adjustCysPosition protSeq cysSeq cysPs a
                elif protSeq.[a-1]='C' then a-1
                elif protSeq.[a]='C' then a
                else adjustCysPosition protSeq cysSeq cysPs a
            
            fillCysItem 0 org protID protSeq cP cysSeq label] )
    |> List.concat
    |> List.toArray
    |> Array.filter (fun i -> i.ProteinSeq.[i.CysPosition]='C')
    
let negativePoolC = 
    itemsC
    |> Array.distinctBy (fun i -> i.ProtID)
    |> Array.map (createNegative itemsC)
    |> Array.concat

let negativePartC =
    let r = new Random()
    Array.sampleWithOutReplacement r (negativePoolC) (itemsC.Length)

let chlamyData =
    [itemsC;negativePartC]
    |> Array.concat
    |> Array.mapi enum
itemsC.Length

/// Read Ara data

let fileDirA = @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\A.fasta"

//reads from file to an array of FastaItems.
let sequencesAra = 
    fileDirA
    |> fromFile (Seq.toArray >> String.fromCharArray)
    |> Seq.toArray
    |> Array.map (fun i -> 
            let h = i.Header |> String.split '|' |> Array.item 0 |> String.subString 0 9//|> String.tryParseIntDefault 0
            {i with Header=h})

let linesAra = 
    [FileIO.readFile @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\A1.txt" |> Seq.toList  |> List.tail;
    FileIO.readFile @"c:\Users\mikha\Work-CSB\Redox-Sensitive Cys\Tables\A2.txt"  |> Seq.toList  |> List.tail;]
    |> List.concat
    |> List.filter (fun i -> i<>"")
    
let itemsAra =
    linesAra 
    |> Seq.toList
    |> List.tail
    |> List.map (fun x -> 
        let i = (x |> String.split '\t')
        let org = i.[1]
        let protID = i.[0]
        printfn "%s" protID
        let protSeq = (sequencesAra |> Array.find (fun i -> i.Header=protID)).Sequence
        let cysPs = i.[2] |> String.replace " ; " ";" |> String.split ';'  |> Array.map (String.replace "\"" "" >> String.tryParseIntDefault 0)
        let label = true
        [for a in cysPs -> fillCysItem 0 org protID protSeq a label] )
    |> List.concat
    |> Array.ofList
    

let negativePoolAra = 
    itemsAra
    |> Array.distinctBy (fun i -> i.ProtID)
    |> Array.map (createNegative itemsAra)
    |> Array.concat

let negativePartAra =
    let r = new Random()
    Array.sampleWithOutReplacement r (negativePoolAra) (itemsAra.Length)

let araData =
    [itemsAra;negativePartAra]
    |> Array.concat
    |> Array.mapi enum
araData.Length

/// the whole data

let data = [redoxData; chlamyData; araData] |> Array.concat |> Array.distinctBy (fun x -> (x.ProteinSeq,x.CysPosition))
data.Length // 8705
itemsC1.Length
araData |> Array.distinctBy (fun x -> (x.ProteinSeq,x.CysPosition)) |> Array.length
data |> Array.filter (fun i -> i.Label=true) |> Array.length
data |> Array.filter (fun i -> i.Label=false) |> Array.length
data |> Array.filter (fun i -> not (i.Organism |> String.contains "Arab") && not (i.Organism |> String.contains "Chlamy")) |> Array.length
data.[0].ProteinSeq


data |> Array.map (fun i -> 
    printfn "seq_Length=%i, item at %i is %c" i.ProteinSeq.Length i.CysPosition i.ProteinSeq.[i.CysPosition]
    i.ProteinSeq.[i.CysPosition])
/////////////// FEATURES