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

open Functions
open Functions.SSN
open Functions.KMeanSwapFunctions
open TestData
open GePhi
open Types
open Auxilliary
open Plots
open MathNet.Numerics
open FunctionsExp
open FunctionsExp
open FSharp.Stats.ML.Unsupervised.HierarchicalClustering


let data6t =
    ChlamyProteome.dataAll
    |> Array.groupBy (fun x -> x.BinL.[0])
    |> Array.map (fun (b,l) -> (b,l.Length, l |> Array.mapi (fun id x -> {x with ID=id})))

let data6t_35 = 
    ChlamyProteome.dataAll
    |> Array.filter (fun x -> x.BinL.[0]="35")
    |> Array.mapi (fun id x -> {x with ID=id})

let matrix_data6t_35 =
    data6t_35 |> SSN.distMatrixWeightedOf SSN.distanceMatrixWeighted None

let data10t =
    ArabiTranscriptome.filteredAra
    |> Array.groupBy (fun x -> x.BinL.[0])
    |> Array.map (fun (b,l) -> (b,l.Length, l |> Array.mapi (fun id x -> {x with ID=id; dataL=x.dataL |> TestData.General.zScoreTransform})))

let data10t_1 = 
    ArabiTranscriptome.filteredAra
    |> Array.filter (fun x -> x.BinL.[0]="1")
    |> Array.mapi (fun id x -> {x with ID=id; dataL=x.dataL |> TestData.General.zScoreTransform})

    
let data10t_13 = 
    ArabiTranscriptome.filteredAra
    |> Array.filter (fun x -> x.BinL.[0]="13")
    |> Array.mapi (fun id x -> {x with ID=id; dataL=x.dataL |> TestData.General.zScoreTransform})


let data10t_26 = 
    ArabiTranscriptome.filteredAra
    |> Array.filter (fun x -> x.BinL.[0]="26")
    |> Array.mapi (fun id x -> {x with ID=id; dataL=x.dataL |> TestData.General.zScoreTransform})

data10t_26.Length

let matrix_data10t_1 =
    data10t_1 |> SSN.distMatrixWeightedOf SSN.distanceMatrixWeighted None

/// 1 Step: apply hierarchical clustering for datapath with ks: [2 .. 2^points]

let pre_k_max_6 = 2*2*2*2*2*2
2.**6.
let pre_k_max_10 = 2*2*2*2*2*2*2*2*2*2

let clusteringSet_6 =
    let basic = 
        data6t_35
        |> Array.map (fun protein -> protein.dataL)
        |> ML.Unsupervised.HierarchicalClustering.generate (weightedEuclidean None) (ML.Unsupervised.HierarchicalClustering.Linker.completeLwLinker)

    [1 .. (min data6t_35.Length pre_k_max_6)] |> List.map (fun k -> 
        (k,basic
        |> ML.Unsupervised.HierarchicalClustering.cutHClust k
        |> List.map (List.map (fun i -> data6t_35.[ML.Unsupervised.HierarchicalClustering.getClusterId i]) >> List.toArray)
        |> List.toArray ))

let clusteringSet_10_fn data =
    let basic = 
        data
        |> Array.map (fun protein -> protein.dataL)
        |> ML.Unsupervised.HierarchicalClustering.generate (weightedEuclidean None) (ML.Unsupervised.HierarchicalClustering.Linker.completeLwLinker)

    [1 .. (min data.Length pre_k_max_10)] |> List.map (fun k -> 
        (k,basic
        |> ML.Unsupervised.HierarchicalClustering.cutHClust k
        |> List.map (List.map (fun i -> data.[ML.Unsupervised.HierarchicalClustering.getClusterId i]) >> List.toArray)
        |> List.toArray ))

let clusteringSet_10_1 = clusteringSet_10_fn data10t_1
let clusteringSet_10_13 = clusteringSet_10_fn data10t_13
let clusteringSet_10_26 = clusteringSet_10_fn data10t_26

/// 2 Step: apply different elbow techniques, plot results, and choose, what suites better.

let plotElbow critTitle (criterion: Item [] [] -> float) (kXclusterings: (int*(Item [] [])) list) =
    let xy = 
        kXclusterings
        |> List.map (fun (k,sch) -> (k, criterion sch))
    let titel = sprintf "Elbow criterion: %s"  critTitle
    Chart.Line (xy, ShowMarkers=true)
    |> Chart.withTitle titel
    |> Chart.Show

let plotElbow' critTitle fn (kXclusterings: (int*(Item [] [])) list) =
    let xy = 
        fn kXclusterings
    let titel = sprintf "Elbow criterion: %s"  critTitle
    Chart.Line (xy, ShowMarkers=true)
    |> Chart.withTitle titel
    |> Chart.Show

plotElbow "Silhouette" (Auxilliary.ClusterCheck.checkSilhouette matrix_data6t_35) clusteringSet_6
plotElbow "Silhouette" (Auxilliary.ClusterCheck.checkSilhouette matrix_data10t_1) clusteringSet_10_1
plotElbow "Silhouette" (Auxilliary.ClusterCheck.checkSilhouette matrix_data10t_1) clusteringSet_10_13
plotElbow "Silhouette" (Auxilliary.ClusterCheck.checkSilhouette matrix_data10t_1) clusteringSet_10_26

let wIndex (data: Item [] []) =
    let centroid (group: Item []) =
        group |> Array.map (fun x -> x.dataL) |> MatrixTopLevelOperators.matrix
    let dist (a: Item) (b: matrix) =
        PreClusterFunctions.pairwiseCorrAverage ([a.dataL] |> MatrixTopLevelOperators.matrix) b
    let w =
        data
        |> Array.map (fun kl -> 
            let c = centroid kl
            kl 
            |> Array.map (fun i ->
                dist i c)
            |> Array.sum) 
        |> Array.sum
    w

plotElbow "Internal W-index" (wIndex) clusteringSet_6
plotElbow "Internal W-index" (wIndex) clusteringSet_10_1
plotElbow "Internal W-index" (wIndex) clusteringSet_10_13
plotElbow "Internal W-index" (wIndex) clusteringSet_10_26

let hartiganIndex (kXclusterings: (int*(Item [] [])) list) =
    let n = kXclusterings.Head |> snd |> Array.concat |> Array.length
    let wk = 
        kXclusterings
        |> List.map (fun (k,sch) -> (wIndex sch))
    [for k in [1 .. (wk.Length-1)] do 
                                        let gamma = float (n-k-1)
                                        yield (k,gamma*(wk.[k-1]-wk.[k])/wk.[k])
    ]

plotElbow' "Hartigan Index" hartiganIndex clusteringSet_6
plotElbow' "Hartigan Index" hartiganIndex clusteringSet_10_1
plotElbow' "Hartigan Index" hartiganIndex clusteringSet_10_13
plotElbow' "Hartigan Index" hartiganIndex clusteringSet_10_26

let sse (data: Item [] []) =
    let centroid (group: Item []) =
        group |> Array.map (fun x -> x.dataL) |> MatrixTopLevelOperators.matrix
    let dist (a: Item) (b: matrix) =
        PreClusterFunctions.pairwiseCorrAverage ([a.dataL] |> MatrixTopLevelOperators.matrix) b
    let w =
        data
        |> Array.map (fun kl -> 
            let c = centroid kl
            kl 
            |> Array.map (fun i ->
                (dist i c)*(dist i c))
            |> Array.average) 
        |> Array.average
    w
 
plotElbow "Sum of Squared Error" (sse) clusteringSet_6
plotElbow "Sum of Squared Error" (sse) clusteringSet_10_1
plotElbow "Sum of Squared Error" (sse) clusteringSet_10_13
plotElbow "Sum of Squared Error" (sse) clusteringSet_10_26

let sseStep (kXclusterings: (int*(Item [] [])) list) =
    let n = kXclusterings.Head |> snd |> Array.concat |> Array.length
    let wk = 
        kXclusterings
        |> List.map (fun (k,sch) -> (sse sch))
    [for k in [1 .. (wk.Length-1)] do
                                        yield (k,(wk.[k-1]-wk.[k]))
    ]

plotElbow' "SSE step" sseStep clusteringSet_6
plotElbow' "SSE step" sseStep clusteringSet_10_1
plotElbow' "SSE step" sseStep clusteringSet_10_13
plotElbow' "SSE step" sseStep clusteringSet_10_26

let dbIndex (data: Item [] []) =
    let centroid (group: Item []) =
        group |> Array.map (fun x -> x.dataL) |> MatrixTopLevelOperators.matrix
    let dist (a: Item) (b: matrix) =
        PreClusterFunctions.pairwiseCorrAverage ([a.dataL] |> MatrixTopLevelOperators.matrix) b
    let s cluster =
        let c = centroid cluster
        cluster
        |> Array.map (fun i ->
                (dist i c)*(dist i c))
        |> Array.average 
        |> sqrt
    let m cl1 cl2 = PreClusterFunctions.pairwiseCorrAverage (centroid cl1) (centroid cl2)
    let sVector = data |> Array.map (fun cl -> s cl)
    let rMatrix = Array2D.init data.Length data.Length (fun i j -> 
                                    if i=j 
                                        then 0. 
                                        else (sVector.[i] + sVector.[j])/(m data.[i] data.[j]))

    rMatrix
    |> Array2D.toJaggedArray
    |> Array.averageBy (fun row -> row |> Array.max)

plotElbow "Davies-Bouldin Index" (dbIndex) clusteringSet_6
plotElbow "Davies-Bouldin Index" (dbIndex) clusteringSet_10_1
plotElbow "Davies-Bouldin Index" (dbIndex) clusteringSet_10_13
plotElbow "Davies-Bouldin Index" (dbIndex) clusteringSet_10_26

let dbIndexStep (kXclusterings: (int*(Item [] [])) list) =
    let n = kXclusterings.Head |> snd |> Array.concat |> Array.length
    let wk = 
        kXclusterings
        |> List.map (fun (k,sch) -> (dbIndex sch))
    [for k in [1 .. (wk.Length-1)] do
                                        yield (k,(wk.[k-1]-wk.[k]))
    ]

plotElbow' "DB-index step" dbIndexStep clusteringSet_6
plotElbow' "DB-index step" dbIndexStep clusteringSet_10_1
plotElbow' "DB-index step" dbIndexStep clusteringSet_10_13
plotElbow' "DB-index step" dbIndexStep clusteringSet_10_26


//// get the point from "monotone" function

let vectorRejection vector projectionLine =
    let pNorm = Vector.norm projectionLine
    let pUnit = Vector.scale (1./pNorm) projectionLine
    let vProj = Vector.scale (Vector.dot vector pUnit) pUnit
    Vector.norm (vector - vProj)

let xy = [1.;-1.] |> MatrixTopLevelOperators.vector

let lineToElbow fn kXclusterings = 
    let xy = kXclusterings |> List.map (fun (k,sch) -> ( float k, fn sch))
    let xmax = xy |> List.map (fst) |> List.max
    let xmin = xy |> List.map (fst) |> List.min
    let ymax = xy |> List.map (snd) |> List.max
    let ymin = xy |> List.map (snd) |> List.min
    xy
    |> List.map (fun (x,y) ->( x,  [(x-xmin)/(xmax-xmin);(y-ymin-ymax)/(ymax-ymin)] |> MatrixTopLevelOperators.vector))

let line = lineToElbow wIndex clusteringSet_10_26 // sse wIndex

let distances = 
    line |> List.map (fun (k,i) -> (k,(i.[0], vectorRejection i xy)))

let optimumFound = distances |> List.maxBy (snd >> snd) |> fst |> int

let optimumGet = (optimumFound/2)+optimumFound

// for (6-timepoints) path 35 - 120 elements: 13 and 19
// for (10-timepoints) path 1 - 150 elements: 37 and 40
// for (10-timepoints) path 13 - 233 elements: 55 and 61
// for (6-timepoints) path 35 - 866 elements: 154 and 210



/// why don't use Gain formula to get the best k???? Try it!