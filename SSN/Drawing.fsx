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
//#load "FunctionsExp.fs"
#load "GePhi.fs"
#load "TestData.fs"
#load "SOM.fs"
#load "Auxilliary.fs"
#load "Plots.fs"

open System
open System.IO
open FSharp.Plotly
open FSharp.Stats

open TestData
open Functions
open Functions.SSN
open Auxilliary
open Plots
//open FunctionsExp

let plain = [1. .. 1. .. 6.]
let recordPoints = [1.;24.;25.;26.;28.;32.]
let recordTimes = [55.;55.;55.;55.;55.;55.] 
let hsLine = [40.;40.]
let hsTime = [0.;24.]
let recLine = [25.;25.]
let recTime = [24.;32.]

[Chart.Line(recTime, recLine, "recovery");
Chart.Line(hsTime, hsLine, "heat acclimation");
Chart.Point(recordPoints, recordTimes, "measured points", Color=colorGray)
]
|> Chart.Combine
|> Chart.withX_AxisStyle("Time, h", (0.,32.5), false, true)
|> Chart.withY_AxisStyle("Temperature, °C", (24.,56.), false, true)
|> Chart.withSize (600.,350.)
|> Chart.Show

let drawLeaves title (tree: Types.Node<string,Types.Item>) =
    tree
    |> Auxilliary.Tree.filterLeaves
    |> (fun x -> x.[0 ..])
    |> Array.mapi (fun c list ->
        let col =
            match c with
            |0 -> colorOrange 
            |1 -> "rgba(0,0,0,1)" 
            |2 -> colorBlue   
            |3 -> colorGreen 
            |4 -> colorYellow 
            |5 -> "rgba(0,180,0,1)"
            |6 -> "rgba(0,180,180,1)"
            |7 -> "rgba(180,0,180,1)"
            |8 -> "rgba(180,180,0,1)"
            |9 -> "rgba(0,0,180,1)"
            |10 -> "rgba(180,0,0,1)"
            |_ -> "rgba(0,0,0,1)"
        list |> Array.map (fun d -> Chart.Line(recordPoints, d.dataL, sprintf "%i" d.ID, Color=col) ) |> Chart.Combine
                )
    |> Chart.Combine
    |> Chart.withTitle title
    |> Chart.withSize (600.,400.)
    |> Chart.Show

let nSet = ChlamyProteome.dataAll.Length

//let data8 = 
//    ChlamyProteome.dataAll
//    |> Array.filter (fun x -> x.BinL.[0]="8")
//    |> Array.mapi (fun id x -> {x with ID=id})

//let data13 = 
//    ChlamyProteome.dataAll
//    |> Array.filter (fun x -> x.BinL.[0]="13")
//    |> Array.mapi (fun id x -> {x with ID=id})

//let data9 = 
//    ChlamyProteome.dataAll
//    |> Array.filter (fun x -> x.BinL.[0]="9")
//    |> Array.mapi (fun id x -> {x with ID=id})
    
let data1 = 
    ChlamyProteome.dataAll
    |> Array.filter (fun x -> x.BinL.[0]="1")
    |> Array.mapi (fun id x -> 
        if id = 40 then
            {x with ID=id; dataL=[|-0.21391304757; 1.621221147; 0.5837022265; -0.2556760858; -0.7524805032; -1.182853737|]}
        else
            {x with ID=id}
        )

//let matrix9 = 
//    data9
//    |> distMatrixWeightedOf distanceMatrixWeighted None
    
//let matrix13 = 
//    data13
//    |> distMatrixWeightedOf distanceMatrixWeighted None

let matrix1 = 
    data1
    |> distMatrixWeightedOf distanceMatrixWeighted None

let treeMMO = Functions.readMM 120 data1
let treeSSN = Functions.applySSNcombi 120 data1

let tree = treeSSN.Children.Item "3"
drawLeaves "Calvin Cycle" tree 

let clustering = [2 .. 7 .. 63] |> List.map (fun k ->  (k,SSN.clusterHier k None (data1 |> Array.toList) |> Map.toArray)) 

GePhi.sendToGephiFromTreeParam treeMMO
GePhi.sendToGephiFromTreeParam treeSSN

let node =
    //treeSSN |> Auxilliary.Tree.findNode ["1";"3";"mix-1-11-2-8"]
    //treeSSN |> Auxilliary.Tree.findNode ["1";"3";"mix-12-4-9"]
    treeSSN |> Auxilliary.Tree.findNode ["1";"3";"mix-13-6"]
    //|> Array.concat

Plots.drawKinetik (node) (recordPoints |> List.toArray)  "path 1.3" |> Chart.Show
drawLeaves "path8" treeSSN

let SSNdxc = Analysis.pointDxC matrix1 (treeSSN |> Tree.filterLeaves)
let MMOdxc = Analysis.SO_points_Fn matrix1 treeMMO

let sortedDxC =
    ChlamyProteome.dataAll
    |> Array.groupBy (fun x -> x.BinL.[0]) 
    |> Array.map (fun (bin,data') ->
        let data = data' |> Array.mapi (fun id x -> {x with ID=id})
        let n = data.Length
        if data.Length > 4 then
            let mmo = data  |> readMM n 
            printfn "mmo for path %s" bin 
            let check = mmo |> Auxilliary.Tree.filterChildrenList |> List.max
            if check < 12 then
                let ssn = applySSNcombi n data
                printfn "draw path %s" bin
                Auxilliary.Analysis.drawDCsinglePath None data mmo ssn
            else
                printfn "not draw path %s" bin
        )


Auxilliary.Analysis.drawDCsinglePath None data1 treeMMO treeSSN

let matrixPath = 
    data1
    |> distMatrixWeightedOf distanceMatrixWeighted (None)

let clustersK =
    [1 .. 5 .. 63] |> List.map (fun k -> 
        data1
        |> Array.map (fun protein -> protein.dataL)
        |> ML.Unsupervised.HierarchicalClustering.generate (weightedEuclidean None) (ML.Unsupervised.HierarchicalClustering.Linker.completeLwLinker)
        |> ML.Unsupervised.HierarchicalClustering.cutHClust k
        |> List.map (List.map (fun i -> data1.[ML.Unsupervised.HierarchicalClustering.getClusterId i]) >> List.toArray)
        |> List.toArray )

let pointSSN = 
    Analysis.pointDxC matrixPath (treeSSN |> Tree.filterLeaves)

let pointsMMO =
    Analysis.SO_points_Fn matrixPath treeMMO

let normalize (maxDissim: float) (input: float*float) =
    let (x,y) = input
    (x/maxDissim,y)
    
let maxDiss = 
    pointsMMO
    |> List.map fst
    |> List.max

let pointSSNnorm = 
    normalize maxDiss pointSSN  

let pointMMOnorm = 
    pointsMMO
    |> List.map (fun i -> normalize maxDiss i) 

let minMMO =
    pointMMOnorm
    |> List.minBy (fun i -> weightedEuclidean None [0.;0.] [fst i; snd i])
    
let pointsHC =
    clustersK
    |> List.map (Analysis.pointDxC matrixPath)
    |> List.map (fun i -> normalize maxDiss i) 

let minHC =
    pointsHC
    |> List.minBy (fun i -> weightedEuclidean None [0.;0.] [fst i; snd i])
    

let drawComparePlot' dataDO (dataSO: (float*float) list) (dataHC: (float*float) list)  =
    
    let doP = 
        Chart.Point ([dataDO], Name = "SSN")
    
    let soP =
        let anno = [0 .. (dataSO.Length-1)] |> List.map string
        Chart.Line (dataSO, Name = "MMO", Labels = anno, ShowMarkers=true)

    let hcP =
        let anno = [0 .. (dataHC.Length-1)] |> List.map string
        Chart.Line (dataHC, Name = "Hierarchical clustering", Labels = anno, ShowMarkers=true)
   
    [soP |> Chart.withLineStyle (Color=colorBlue);
    doP |> Chart.withLineStyle (Color=colorOrange);
    hcP |> Chart.withLineStyle (Dash = StyleParam.DrawingStyle.Dash, Color=colorGray);
    Chart.Line [(0.0, 0.0); pointSSNnorm] |> Chart.withLineStyle (Color=colorOrange);
    Chart.Line [(0.0, 0.0); minMMO] |> Chart.withLineStyle (Color=colorBlue);
    Chart.Line [(0.0, 0.0); minHC] |> Chart.withLineStyle (Color=colorGray);
    ]
    |> Chart.Combine
    |> Chart.withLegend (false)
    |> Chart.withX_AxisStyle ("Purity, norm", Showgrid=false)
    |> Chart.withY_AxisStyle ("Purity, norm", Showgrid=false)
    |> Chart.withSize (500.,500.)
    |> Chart.Show

drawComparePlot' pointSSNnorm pointMMOnorm pointsHC

let mmoP = sqrt((0.2857143*0.2857143)+(0.3423718698*0.3423718698))
let ssnP = sqrt((0.1133837428*0.1133837428)+(0.5079365079*0.5079365079))
let hccP = sqrt((0.253968254*0.253968254)+(0.2210243002*0.2210243002))