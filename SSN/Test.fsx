#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\netstandard.dll"

#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\BioFSharp.dll"
#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\BioFSharp.IO.dll"
#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharpAux.dll"
#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharpAux.IO.dll"
#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\Newtonsoft.Json.dll"
#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharp.Stats.dll"
#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharp.Plotly.dll"
#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharpGephiStreamer.dll"

#r @"c:\Users\mikha\source\repos\mathnet-numerics\src\Numerics\bin\Debug\netstandard2.0\MathNet.Numerics.dll"

#load "Types.fs"
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

open Functions
open FunctionsExp
open TestData
open GePhi
open Types
open Auxilliary
open FunctionsExp
open FunctionsExp
open FunctionsExp
open FSharp.Stats.ML

/// comparison of pairwise matrices (each is created as 0 if the pair is in the same leaf, 1 if not). O(n*n)
let treeComparison (treeA: Node<string,Item>) (treeB: Node<string, Item>) : int =
    let m tree =
        tree 
        |> Tree.filterLeaves |> Array.concat |> Array.map (fun x -> x.ID,x.BinL) |> Array.sortBy fst
        |> Array.pairwise |> Array.map (fun ((id1,bin1),(id2,bin2)) -> if bin1=bin2 then 1 else 0)
    let lA = m treeA 
    let lB = m treeB 
    Array.zip lA lB |> Array.sumBy (fun ((a),(b)) -> if a=b then 0 else 1)

//let dataInfo =
//    (ArabiTranscriptome.filteredAra) //ChlamyProteome.dataAll
//    |> Array.groupBy (fun x -> x.BinL.[0])
//    |> Array.map (fun (bin,l) -> (bin, l|> Array.length))
//    |> Array.sortBy (fst >> int)

//let nSet = (ArabiTranscriptome.filteredAra).Length

//let childrenInfo =
//    (ArabiTranscriptome.filteredAra)
//    |> Array.groupBy (fun x -> x.BinL.[0])
//    |> Array.filter (fun (_,l) -> l.Length>2)
//    |> Array.map (fun (bin,items') -> 
//        let items =
//            items'
//            |> Array.mapi (fun id it -> {it with ID=id; dataL = it.dataL |> General.zScoreTransform})
//        let childrenSize = (readMM items nSet |> Tree.filterChildrenList |> List.toArray |> Array.sortDescending)
        
//        (bin, childrenSize)
        
//        )

//childrenInfo |> Array.sortBy (fun (_,c) -> c .[0])

let pathsSorted =
    ChlamyProteome.dataAll
    |> Array.groupBy (fun x -> x.BinL.[0])
    |> Array.filter (fun (bin,l) -> l.Length > 2)
    |> Array.sortBy (fun (bin,l) -> l.Length )
    |> Array.map (fst)

let data16 = 
    ChlamyProteome.dataAll
    |> Array.filter (fun x -> x.BinL.[0]="16")
    |> Array.mapi (fun id x -> {x with ID=id})

let data29 = 
    ChlamyProteome.dataAll
    |> Array.filter (fun x -> x.BinL.[0]="29")
    |> Array.mapi (fun id x -> {x with ID=id})

let data35 = 
    ChlamyProteome.dataAll
    |> Array.filter (fun x -> x.BinL.[0]="35")
    |> Array.mapi (fun id x -> {x with ID=id})


let data35toPre = 
    ChlamyProteome.dataAll
    |> Array.filter (fun x -> x.BinL.[0]="35" && x.BinL.[1]="2")
    |> Array.mapi (fun id x -> {x with ID=id})

let data35toPreMap = 
    data35toPre
    |> Array.map (fun x -> (string x.ID, [|x|]))
    |> Map.ofArray

let matrix = 
    data35toPre
    |> SSN.distMatrixWeightedOf SSN.distanceMatrixWeighted None

let newClusters = PreClusterFunctions.clusterHierGroups 30 None data35toPreMap

let gainFn = SSN.getStepGainNodeSetnR 300
let preG = 
    newClusters 
    |> Map.toArray 
    |> Array.sumBy (fun (s,items) -> SSN.getStepGainFn gainFn (items |> SSN.groupIDFn) (data35toPre |> SSN.groupIDFn) data35toPre.Length matrix)
printfn "pre-Cluster G = %f" preG

let mmPQ = FunctionsExp.readMMpq data29 120 |> Tree.filterChildrenList
let mmNoPQ = FunctionsExp.readMM data35 120 |> Tree.filterChildrenList
//mmPQ = mmNoPQ

let ssnNoPQ = FunctionsExp.applySSN data29 300                      // |> ignore // 9.7 G=120.8      
let ssnPreClus = FunctionsExp.applySSN_Reduction 30 data29 300         // |> ignore // 8.7 G=120.8   
let ssnPreClusR = FunctionsExp.applySSN_Reduction_Random 30 data29 300 // |> ignore // 7.4 G=120.8   

let ssnOld35 = FunctionsExp.applySSN data35 300                         // |> ignore // 18:23 G=84.07421, compare it to G=65 for pre_k=50
let ssnPreClus35 = FunctionsExp.applySSN_Reduction 30 data35 300           // |> ignore // 00:22 G=78.3 with pre_k=30 with random swap G=84.6
let ssnPreClusR35 = FunctionsExp.applySSN_Reduction_Random 30 data35 300   // |> ignore // 00:20 G=81.1 ; 78.6 ; 76.9 ; 77.2 ; 86.7!!! 93!!!
let ssnPreClusH35 = FunctionsExp.applySSN_Reduction_Hier 30 data35 300   // |> ignore // 00:20 G=81.1 ; 78.6 ; 76.9 ; 77.2 ; 86.7!!! 93!!!

data35 |> Array.map (fun i -> i.dataL |> Array.sum)

let f =
    [for a in [1 .. 50] ->    
        let x = FunctionsExp.applySSN_Reduction_Random 30 data35 300
        printfn "G=%f" x.GroupGain
        x
        ]

let best = f |> List.maxBy (fun i -> i.GroupGain)
best.GroupGain

treeComparison ssnPreClus35 ssnPreClusR35

ssnNoPQ |> Tree.filterChildrenList
ssnPreClus |> Tree.filterChildrenList


let ssnReducedComplexity = FunctionsExp.applySSN_Reduction 30 data35 300               

ssnReducedComplexity |> Tree.filterChildrenList

ssnNoPQ.ConfGain
 //////////

GePhi.sendToGephiFromTreeParam best

let node =
    best
    |> Auxilliary.Tree.findNode ["35";"2";"||p51|p83|p47|p50|p17|p84|p66||p49|p97|p1|p10||p22|p81|p120|p102|p58|p8||p0|p74|p39|p96||p20|p35|p93||p107|p75||p100|p37|p15|p29|p34||p61|p64|p28|p95|p48|p57|p63|p108|p82"]

Plots.drawKinetik (node) [|1.;24.;25.;26.;28.;32.|]  "subbin 35.2.mix big" |> Chart.Show

let allPathsWriteNormal listN items =
        //let pathFile = @"c:\_n_mikhaylenko\Code_FSharp\Projects\SSN\results\check\"
        //let pathFile = @"D:\Nathan\results\tablesAraArraysOnly\"
        let pathFile = @"c:\Users\mikha\Work-CSB\Projects\SSN\results\checkNoPQ\"

        listN
        |> Array.mapi (fun i n ->
                    
                        let itemsOfN =
                            items
                            |> Array.filter (fun i -> i.BinL.Length>0 && i.BinL.[0] = n)
                            |> Array.mapi (fun index i -> {i with ID=index; dataL = i.dataL |> General.zScoreTransform })

                        let treeMMofN = FunctionsExp.readMM itemsOfN items.Length
                        let stopwatch = new System.Diagnostics.Stopwatch()
                        stopwatch.Start()
                        let treeSSNofN = FunctionsExp.applySSN itemsOfN items.Length                //// change here which SSN to use
                        let time = sprintf "SSN tree calculated in %f s since start" (stopwatch.Elapsed.TotalSeconds)
                        stopwatch.Stop()
                        let dc = Auxilliary.Analysis.getDCmeasure None itemsOfN treeMMofN treeSSNofN
                        let fileName = n
                        let title = sprintf "bin: %s, items: %i" n itemsOfN.Length
                        let dcText = sprintf "DxC measure: %f" dc
                        let header = "Bin\tItem"
                        let content = 
                            treeSSNofN 
                            |> Auxilliary.Tree.mapTreeToBinsMembers 
                            |> Array.toList 
                            |> List.map (fun (bin,x) -> sprintf "%s\t%A" (String.Join(".", bin)) x.ProteinL)
                        //treeArray.[i] <- fastTreeOfN
                        File.WriteAllLines((sprintf "%s%s.txt" pathFile fileName), title :: time :: dcText :: header :: content)
                        )

let allPathsWritePQ listN items =
        //let pathFile = @"c:\_n_mikhaylenko\Code_FSharp\Projects\SSN\results\check\"
        //let pathFile = @"D:\Nathan\results\tablesAraArraysOnly\"
        let pathFile = @"c:\Users\mikha\Work-CSB\Projects\SSN\results\checkPQnoCopy\"

        listN
        |> Array.mapi (fun i n ->
                    
                        let itemsOfN =
                            items
                            |> Array.filter (fun i -> i.BinL.Length>0 && i.BinL.[0] = n)
                            |> Array.mapi (fun index i -> {i with ID=index; dataL = i.dataL |> General.zScoreTransform })

                        let treeMMofN = FunctionsExp.readMMpq itemsOfN items.Length
                        let stopwatch = new System.Diagnostics.Stopwatch()
                        stopwatch.Start()
                        let treeSSNofN = FunctionsExp.applySSNpq itemsOfN items.Length                //// change here which SSN to use
                        let time = sprintf "SSN tree calculated in %f s since start" (stopwatch.Elapsed.TotalSeconds)
                        stopwatch.Stop()
                        let dc = Auxilliary.Analysis.getDCmeasure None itemsOfN treeMMofN treeSSNofN
                        let fileName = n
                        let title = sprintf "bin: %s, items: %i" n itemsOfN.Length
                        let dcText = sprintf "DxC measure: %f" dc
                        let header = "Bin\tItem"
                        let content = 
                            treeSSNofN 
                            |> Auxilliary.Tree.mapTreeToBinsMembers 
                            |> Array.toList 
                            |> List.map (fun (bin,x) -> sprintf "%s\t%A" (String.Join(".", bin)) x.ProteinL)
                        //treeArray.[i] <- fastTreeOfN
                        File.WriteAllLines((sprintf "%s%s.txt" pathFile fileName), title :: time :: dcText :: header :: content)
                        )

allPathsWriteNormal pathsSorted ChlamyProteome.dataAll
allPathsWritePQ [|"35";"29"|] ChlamyProteome.dataAll

[1.;1.5;1.2;2.;4.;8.;2.1] |> List.median

//// Centroid test
//let points = [3.,5.;7.,2.;-5.,1.;-6.,-2.;-3.,-5.]

let points = [3.,7.;3.,3.;5.,5.;5.,3.;7.,3.]

let centroidMedian = (points |> List.map (fst) |> List.median |> fun (a,b) -> (a+b)/2.),(points |> List.map (snd) |> List.median |> fun (a,b) -> (a+b)/2.)
let centroidAverage = (points |> List.map (fst) |> List.average),(points |> List.map (snd) |> List.average)

let dist (x0,y0) =
    points |> List.sumBy (fun (x,y) -> DistanceMetrics.cityblock [x0;y0] [x;y])
let minimumA = dist centroidAverage
let minimumM = dist centroidMedian

let set = 
    [for x in [-10.0 .. 0.05 .. 10.] do 
        for y in [-10. .. 0.05 .. 10.] do 
            yield x,y,(dist (x,y))]
            //printfn "for (%f,%f), dist = %f" x y (dist (x,y))]

let lower = set |> List.filter (fun (x,y,d) -> d<(min minimumA minimumM)) |> List.map (fun (x,y,_) -> (x,y))
let bestC = set |> List.minBy (fun (x,y,d) -> d) |> (fun (x,y,_) -> (x,y))
dist bestC

[
Chart.Scatter(points, StyleParam.Mode.Markers, "Dataset");
Chart.Scatter(lower, StyleParam.Mode.Markers, "Better centroids");
Chart.Scatter([bestC], StyleParam.Mode.Markers, "best point, SUM(d)=28.79");
Chart.Scatter([centroidMedian], StyleParam.Mode.Markers, "centroid Median, SUM(d)=29.84");
Chart.Scatter([centroidAverage], StyleParam.Mode.Markers, "centroid Average, SUM(d)=30.2");
]
|> Chart.Combine
|> Chart.withSize (800.,600.)
|> Chart.Show