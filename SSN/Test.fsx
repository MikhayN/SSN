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

#load "Types.fs"
#load "PQ.fs"
#load "Functions.fs"
#load "FunctionsExp.fs"
#load "GePhi.fs"
#load "TestData.fs"
#load "SOM.fs"
#load "Auxilliary.fs"
#load "Plots.fs"

#time

open System
open System.IO
open FSharp.Stats
open FSharp.Plotly
open BioFSharp.IO
open FSharpAux.IO
open FSharpAux

open Functions
open TestData
open GePhi
open Types
open Auxilliary
open FSharp.Stats.ML

let readStringBin str =
    str 
    |> String.replace "\"" "" |> String.replace "[|" "" |> String.replace "|]" "" |> String.replace " " ""
    |> String.split ';'

let path1_182_items =
    FileIO.readFile (@"c:\Users\mikha\Downloads\path1_182.txt")
    |> Seq.toArray
    |> Array.map (String.split '\t')
    |> Array.map (fun i -> {ID= String.tryParseIntDefault -1 i.[0]; ProteinL=[||]; OriginalBin = readStringBin i.[1]; BinL = readStringBin i.[2] ; dataL=[||] })
    |> Array.sortBy (fun x -> x.ID)

let dataPOI = 
    ArabiProteome.itemsWithMapManFound
    |> Array.filter (fun x -> x.BinL.[0]="1")
    |> Array.mapi (fun id x -> {x with ID=id})

let dataInTree =
    dataPOI
    |> Array.map (fun i -> 
        let findBin = path1_182_items.[i.ID].BinL
        {i with BinL=findBin})

let tree = readMM_raw dataInTree.Length dataInTree

dataInTree.[25] // "ATCG00490" // ID=25, Rubisco
dataInTree.[90] // "AT2G39730" // ID=90, Rubisco activase

sendToGephiFromTreeParam tree

Plots.drawKinetikRangeStack [|0. .. 7.|] "" (tree |> Tree.findNode ["1";"3"] |> Tree.filterLeaves) |> Chart.Show