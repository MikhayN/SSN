﻿#r @"..\lib\MathNet.Numerics.dll"
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

//open Functions
open Functions
open Functions.SSN
open TestData
open GePhi
open Types
open PQ
open Auxilliary
open Plots
//open FunctionsExp

#load "ElbowCriteria.fs"
open ElbowCriteria

#time

/// cluster elements with Hierarchical clustering and k found from Elbow criterion  
let groupClust data =
    data
    |> GetK.optimumFound ElbowCrit.wIndex
    |> snd

/// get true classes for reference
let trueClasses data =
    data
    |> Array.groupBy (fun i -> i.OriginalBin)
    |> Array.map (snd)

let getSST data =
    data
    |> applySST_walk data.Length
    |> Tree.filterLeaves

//// Parameters

/// purity of cluster set, with reference of OriginalBin as true classes
let purity (clusterSet: Item [] []) n =
    clusterSet
    |> Array.sumBy (fun iList -> 
        iList 
        |> Array.groupBy (fun i -> i.OriginalBin) 
        |> Array.maxBy (fun (c,n_ck) -> n_ck.Length) 
        |> snd 
        |> Array.length)
    |> float
    |> fun s -> s/(float n)

/// complexity of cluster set
let complexity (clusterSet: Item [] []) n =
    (float clusterSet.Length)/(float n)

/// dissimilarity as a n averange cluster multidimensional sd
let dissimilarity (clusterSet: Item [] []) n =
    if n>6 then
        let data cluster =
            cluster
            |> Array.map (fun i -> i.dataL)

        let meanV (clusterData: float [] []) =
            clusterData
            |> Array.transpose
            |> Array.map (Array.average)
    
        clusterSet
        |> Array.averageBy (fun cluster ->
            let d = cluster |> data
            let m = meanV d
            sqrt((d |> Array.sumBy (fun x -> square (weightedEuclidean None x m) ))/(float(n-6)))
            )

    else 
        printfn "dataset size is too small (must be bigger than kinetic vector length)"
        0.

/// similarity between cluster means and classes means, expressed in heatmap plot
let similarity (clusterSet: Item [] []) (classesSet: Item [] []) =
    let meanV (clusterData: float [] []) =
            clusterData
            |> Array.transpose
            |> Array.map (Array.average)
    
    let data =
        clusterSet |> Array.map (fun ik ->
            classesSet |> Array.map (fun jc ->
                let i = ik |> Array.map (fun x -> x.dataL) |> meanV
                let j = jc |> Array.map (fun x -> x.dataL) |> meanV
                Correlation.Seq.pearson i j
            ))

    Chart.Heatmap (data, [1 .. clusterSet.Length], [1 .. classesSet.Length], "Correlation between clusters and classes", true)
    

////// TEST


let pathsSortedChildren =
    ChlamyProteome.dataAll  //ChlamyProteome.dataAll //ArabiTranscriptome.itemsWithMapManIdentified //ArabiProteome.itemsWithMapManFound 
    |> Array.groupBy (fun x -> x.BinL.[0])
    |> Array.filter (fun (bin,l) -> l.Length > 2)
    |> Array.map (fun (bin,l) -> 
        let data = l |> Array.mapi (fun id x -> {x with ID=id})
        let tree = readMM data.Length data 
        let childrenN = Tree.filterChildrenList tree
        (bin, childrenN |> List.max, tree.GroupGain) )
    |> Array.sortBy (fun (_,cn,_) -> cn)
    |> Array.map (fun (x,_,_) -> x)


let dataPOI = 
    ChlamyProteome.dataAll
    |> Array.filter (fun x -> x.BinL.[0]="1")
    |> Array.sortBy (fun x -> x.ProteinL)
    |> Array.mapi (fun id x -> {x with ID=id})

dataPOI

dataPOI.Length

let classes = trueClasses dataPOI

classes.Length

let rootClusters =
    dataPOI
    |> groupClust

rootClusters |> Array.map Array.length 

let leavesClusters =
    classes
    |> Array.map groupClust
    |> Array.concat

leavesClusters.Length

let clustersSST = getSST dataPOI

clustersSST.Length

purity classes (dataPOI.Length) // 1.
purity rootClusters (dataPOI.Length) // 0.46
purity leavesClusters (dataPOI.Length) // 1.
purity clustersSST (dataPOI.Length) // 0.54

complexity rootClusters (dataPOI.Length) // 0.19
complexity classes (dataPOI.Length) // 0.38
complexity leavesClusters (dataPOI.Length) // 0.62
complexity clustersSST (dataPOI.Length) // 0.48

dissimilarity rootClusters (dataPOI.Length) // 0.181
dissimilarity classes (dataPOI.Length) // 0.185
dissimilarity leavesClusters (dataPOI.Length) // 0.04
dissimilarity clustersSST (dataPOI.Length) // 0.07

////// find a way to rearrange the cells in heatmap plot
similarity rootClusters classes |> Chart.Show //
similarity rootClusters rootClusters |> Chart.Show //

similarity leavesClusters classes |> Chart.Show //
similarity leavesClusters leavesClusters |> Chart.Show //


similarity clustersSST classes |> Chart.Show //


let dataSet =
    ChlamyProteome.dataAll
    |> Array.groupBy (fun x -> x.BinL.[0])
    |> Array.map (snd >> Array.mapi (fun id x -> {x with ID=id}))
    |> Array.filter (fun il -> il.Length>6 && il.[0].OriginalBin.[0]<>"35")

dataSet |> Array.iter (fun i -> printfn "%s" i.[0].OriginalBin.[0])

let writeTreeInFile len data = async {
    let tree = applySST_walk len data
    let lines =
        tree 
        |> Tree.filterLeaves 
        |> Array.concat 
        |> Array.map (fun i -> sprintf "%i\t%s\t%s" i.ID (String.Join(";", i.ProteinL)) (String.Join(".", i.BinL)) )
    let x = (File.AppendAllLines((sprintf "%sresults\\SST_%s.txt" General.pathToData data.[0].BinL.[0]), lines))
    return tree
    }

let trees =
    dataSet
    |> Array.map (fun x -> writeTreeInFile x.Length x)
    |> Async.Parallel
    |> Async.RunSynchronously
