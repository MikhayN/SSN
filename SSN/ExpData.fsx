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

open Functions
open Functions.SSN
open Functions.KMeanSwapFunctions
open TestData
open GePhi
open Types
open Auxilliary
open Plots
open MathNet.Numerics

//let fastTreePure = applySSNcombi pureSignal 40

let data6t =
    ChlamyProteome.dataAll
    |> Array.groupBy (fun x -> x.BinL.[0])
    |> Array.map (fun (b,l) -> (b,l.Length, l |> Array.mapi (fun id x -> {x with ID=id})))

let data10t =
    ArabiTranscriptome.filteredAra
    |> Array.groupBy (fun x -> x.BinL.[0])
    |> Array.map (fun (b,l) -> (b,l.Length, l |> Array.mapi (fun id x -> {x with ID=id; dataL=x.dataL |> TestData.General.zScoreTransform})))

let pattern (values: 'a list) (power: int) : 'a list list =
    let rec loop p acc =
        [for a in values do
            if p=1
                then yield (a::acc)
                else yield! loop (p-1) (a::acc)]
    loop power []

pattern [0;1] 2

let filterPatterns l data =
    let patterns = pattern [0;1] l
    let binData = 
        data
        |> Array.map (fun x -> x.dataL |> Array.map (fun i -> if i>0. then 1 else 0)) 
        |> Array.distinct
    printfn "patterns: %i" patterns.Length
    printfn "in data: %i" binData.Length
    (float binData.Length) / ( float patterns.Length)

data6t |> Array.map (fun (b,l,x) -> (b,l,filterPatterns 6 x)) |> Array.sortByDescending (fun (b,l,p) -> p)
data10t |> Array.map (fun (b,l,x) -> (b,l,filterPatterns 10 x)) |> Array.sortByDescending (fun (b,l,p) -> p)
