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

open FSharp.Stats.Correlation

open FSharp.Stats.ML.DistanceMetrics

open FSharp.Stats.ML
open FSharp.Stats.ML.Unsupervised

open MathNet.Numerics.LinearAlgebra.Double
open MathNet.Numerics.LinearAlgebra.Factorization
open MathNet.Numerics.LinearAlgebra

//open FSharp.Stats.Algebra
//open FSharp.Stats.Correlation

#time

//let singleton =
//    [[1.;0.;0.;0.;-2.1]]
//    //|> MatrixTopLevelOperators.matrix 

//let mTest =
//    [|[|1.;0.;0.;0.;2.|];
//    [|0.;0.;3.;0.;0.|];
//    [|0.;0.;0.;0.;0.|];
//    [|0.;2.;0.;0.;0.|]|]
//    |> DenseMatrix.OfRowArrays
    
//mTest.Svd().VT.Row 0

//let mTestClose =
//    [[1.;0.;0.;0.;2.1];
//    [0.;0.;0.;0.;0.];
//    [0.;0.;3.;0.;0.];
//    [0.;2.1;0.;0.;0.]]
//    |> MatrixTopLevelOperators.matrix 

//let sumM =
//    [[1.;0.;0.;0.;2.];
//    [0.;0.;3.;0.;0.];
//    [0.;0.;0.;0.;0.];
//    [0.;2.;0.;0.;0.];
//    [1.;0.;0.;0.;2.1];
//    [0.;0.;3.;0.;0.];
//    [0.;0.;0.;0.;0.];
//    [0.;2.1;0.;0.;0.]]
//    |> MatrixTopLevelOperators.matrix 


//let mTestFar =
//    [[1.;0.;0.;0.;-2.1];
//    [0.;0.;3.;0.;-3.];
//    [10.;0.;0.;0.;0.];
//    [3.;0.;0.;-3.;0.]]
//    //|> MatrixTopLevelOperators.matrix 

//LinearAlgebra.SVD mTest
//LinearAlgebra.SVD mTestClose
//LinearAlgebra.SVD mTestFar
//LinearAlgebra.SVD sumM
//LinearAlgebra.SVD singleton.Transpose


/////// Synthetic data

//let f = (getStepGainNodeSetnR 20)

//let generateSyn pattern noiseSigma n =
//    let r = MathNet.Numerics.Distributions.Normal(0., noiseSigma)
//    [|1 .. n|] 
//    |> Array.mapFold (fun prev i -> (pattern*log(prev*(float i))),(prev+1.)) 1.
//    |> fst 
//    |> General.zScoreTransform
//    |> Array.map (fun x -> x + (r.Sample()))
//    |> General.zScoreTransform

//let generateSynDiff pattern noiseSigma n =
//    let r = MathNet.Numerics.Distributions.Normal(0., noiseSigma)
//    [|1 .. n|] 
//    |> Array.mapFold (fun prev i -> (pattern*log((abs pattern)*prev*(float i))),(prev+1.)) 1. //
//    |> fst 
//    |> zScoreTransform
//    |> Array.map (fun x -> x + (r.Sample()))
    
//let synData1Level pattern noise x =
//    let p = 
//        if pattern > 0. then 1
//        else 2    
//    Array.init x (fun i -> 
//        {ID = i;
//        ProteinL = [|sprintf "name-%i" p|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL = generateSyn pattern noise 6};)

//let synData2Levels pattern noise x =
//    let p = 
//        if pattern > 0. then 1
//        else 2  
//    Array.init x (fun i ->
//        {ID = i;
//        ProteinL = [|sprintf "name-%i" p|];
//        OriginalBin = [|"root";sprintf "%f" pattern|];
//        BinL = [|"root";sprintf "%f" pattern|];
//        dataL = generateSyn pattern noise 6};)

//let pureSignal = 
//    Array.append (synData1Level -1. 0. 2) (synData1Level 1. 0. 2) 
//    |> Array.mapi (fun i p -> {p with ID=i})  

//let noisySignal noise = 
//    Array.append (synData1Level 1. noise 5) (synData1Level -1. noise 5) 
//    |> Array.mapi (fun i p -> {p with ID=i;ProteinL=[|sprintf "%s-%i" p.ProteinL.[0] i|]})

//let noisySignal2 = 
//    Array.append (synData2Levels 1. 0.1 10) (synData2Levels -1. 0.1 10) 
//    |> Array.mapi (fun i p -> {p with ID=i;ProteinL=[|sprintf "%s-%i" p.ProteinL.[0] i|]})

//let noisySignal2ins = 
//    noisySignal2 
//    |> Array.mapi (fun id p -> if id=1 || id=3 || id=5 || id=6 then {p with OriginalBin = [|"root"; "-1"|]; BinL = [|"root"; "-1"|];} else p) 
//    |> Array.mapi (fun id p -> if id=11 || id=13 || id=15 || id=16 then {p with OriginalBin = [|"root"; "1"|]; BinL = [|"root"; "1"|];} else p)

//// create plot with Synthetic data

////drawKinetik (synData1Level 1. 0.1 2) [|1.; 2.; 3.; 4.; 5.; 6.|] "noised signal, up" |> Chart.withLineStyle(Color="#1f77b4");
////drawKinetik (synData1Level -1. 0.1 2) [|1.; 2.; 3.; 4.; 5.; 6.|] "noised signal, down" |> Chart.withLineStyle(Color="#ff7f0e");
////drawKinetik (synData1Level 1. 0.3 2) [|1.; 2.; 3.; 4.; 5.; 6.|] "noised signal, up" |> Chart.withLineStyle(Color="#1f77b4");
////drawKinetik (synData1Level -1. 0.3 2) [|1.; 2.; 3.; 4.; 5.; 6.|] "noised signal, down" |> Chart.withLineStyle(Color="#ff7f0e");
//[drawKinetik (synData1Level 1. 1.1 2) [|1.; 2.; 3.; 4.; 5.; 6.|] "noised signal, up" |> Chart.withLineStyle(Color="#1f77b4");
//drawKinetik (synData1Level -1. 1.1 2) [|1.; 2.; 3.; 4.; 5.; 6.|] "noised signal, down" |> Chart.withLineStyle(Color="#ff7f0e");
//Chart.Spline([1.; 2.; 3.; 4.; 5.; 6.],(generateSyn 1. 0. 6),Name="pattern, up",Color="#144c72", Dash=StyleParam.DrawingStyle.Dash, Width=3.);
//Chart.Spline([1.; 2.; 3.; 4.; 5.; 6.],(generateSyn -1. 0. 6),Name="pattern, down",Color="#ac560c", Dash=StyleParam.DrawingStyle.Dash, Width=3.);
//]
//|> Chart.Combine
//|> Chart.withX_AxisStyle ("Time", Showgrid = false, Showline = true)
//|> Chart.withY_AxisStyle ("Abundance", (-2.,2.), Showgrid= false, Showline = true)
//|> Chart.withSize (650., 300.)
//|> Chart.Show


//drawKinetik (noisySignal 0.3) [|1.; 2.; 3.; 4.; 5.; 6.|] "noise level: 0.3"

//let fastTreePure = applySSNcombi pureSignal 40
//let fastTreeNoise = applySSNcombi (noisySignal 0.5) 60 
//let fastTreeNoise2 = applySSN noisySignal2 40 
//let fastTreeNoise2ins = applySSN noisySignal2ins 40  

//let ss = pureSignal |> Array.map (fun p -> p.dataL) |> Array.concat

//let empDist = FSharp.Stats.Distributions.Empirical.create 1. (ss)
//Chart.Spline (empDist |> Map.toSeq) |> Chart.Show
//let sigmaSignal = FSharp.Stats.Distributions.Empirical.var empDist

/////// 1 setup: noise robustness

///// error as a wrong mixed elements
//let errorFnMix fastTreeNoise = 
//    fastTreeNoise.Children 
//    |> Map.toArray 
//    |> Array.sortByDescending (fun (k,node) -> node.Member.Length)
//    |> Array.mapi (fun id (k,node) ->
//        let sizeList = 
//            node.Member 
//            |> Array.groupBy (fun p -> String.subString 5 1 p.ProteinL.[0]) 
//            |> Array.sortByDescending (fun (name,l) -> l.Length) 
//            |> Array.map (snd >> Array.length)
//        if sizeList.Length>1 
//            then sizeList.[1]
//            else 0
//        )
//    |> Array.sum

///// error as an unnecessary group
//let errorFnAdd fastTreeNoise = 
//    (fastTreeNoise.Children 
//    |> Map.count) - 2

//let errorFn fastTreeNoise = 
//    fastTreeNoise.Children 
//    |> Map.toArray 
//    |> Array.sortByDescending (fun (k,node) -> node.Member.Length)
//    |> Array.mapi (fun id (k,node) ->
//        let sizeList = 
//            node.Member 
//            |> Array.groupBy (fun p -> String.subString 5 1 p.ProteinL.[0]) 
//            |> Array.sortByDescending (fun (name,l) -> l.Length) 
//            |> Array.map (snd >> Array.length)
//        if id<2 then 
//            if sizeList.Length>1 
//                then sizeList.[1]
//                else 0
//        else 
//            sizeList 
//            |> Array.sum)
//    |> Array.sum

//let checkStat =
//    [0. .. 0.1 .. 1.5]
//    |> List.map (fun noiseS ->
//        printfn "noise level: %f" noiseS
//        let averError =
//            [for x in 1 .. 10 ->
//                printfn "iteration %i" x
//                let noisySignal = 
//                    Array.append (synData1Level 1. noiseS 5) (synData1Level -1. noiseS 5)
//                    |> Array.mapi (fun i p -> {p with ID=i})
//                let fastTreeNoise = applySSNcombi noisySignal 40  
//                let error = errorFn fastTreeNoise
//                let mixE = errorFnMix fastTreeNoise
//                let mixA = errorFnAdd fastTreeNoise
//                printfn "error: %i + %i = %i" mixE mixA error
//                (float mixE,float mixA)
//            ]
//            //|> List.average
//        (noiseS,averError)
//    )

//let xy = checkStat |> List.map (fun (x,yyy) -> yyy |> List.map (fun y -> (string x,y))) |> List.concat |> List.toSeq
//let line = checkStat |> List.map (fun (x,yyy) -> (x, yyy |> List.averageBy (fun (e1,e2) -> e1+e2)))
//let (x,y) = xy |> Seq.unzip
//let (e1,e2) = y |> Seq.unzip

//[Chart.BoxPlot(x,e2); Chart.BoxPlot(x,e1); Chart.Line line]
//|> Chart.Combine
//|> Chart.withX_AxisStyle("noise level", Showgrid=false, Showline=false)
//|> Chart.withY_AxisStyle("errors", (0.,6.1), Showgrid=true)
//|> Chart.Show

////// 2. setup: depth vs noise level

//let template = [|
//    {ID = 0;
//    ProteinL = [|"Cre08.g372950.t1.1"|];
//    OriginalBin = [|"16"; "1"; "1"; "7"|];
//    BinL = [|"16"; "1"; "1"; "7"|];
//    dataL = [||];};
//    {ID = 1;
//    ProteinL = [|"Cre12.g509650.t1.1"|];
//    OriginalBin = [|"16"; "2"|];
//    BinL = [|"16"; "2"|];
//    dataL = [||];};
//    {ID = 2;
//    ProteinL = [|"Cre03.g207800.t1.1"|];
//    OriginalBin = [|"16"; "1"; "4"|];
//    BinL = [|"16"; "1"; "4"|];
//    dataL = [||];};
//    {ID = 3;
//    ProteinL = [|"Cre01.g050950.t1.1"|];
//    OriginalBin = [|"16"; "1"; "1"|];
//    BinL = [|"16"; "1"; "1"|];
//    dataL = [||];};
//    {ID = 4;
//    ProteinL = [|"Cre07.g356350.t1.1"|];
//    OriginalBin = [|"16"; "1"; "1"; "1"|];
//    BinL = [|"16"; "1"; "1"; "1"|];
//    dataL = [||];}|]

//let fillProtein bin pattern noise =
//    let name = if pattern>0. then [|"uppp"|] else [|"down"|]
//    {ID = 0;
//    ProteinL = name;
//    OriginalBin = bin;
//    BinL = bin;
//    dataL = generateSyn pattern noise 6;}


//let synthDeepData noise id= 
//    template
//    |> Array.mapi (fun index i -> 
//        if index=id then
//            {i with dataL=generateSyn -1. noise 6}
//        else
//            {i with dataL=generateSyn 1. noise 6} )

//let checkStat2 =
//    [0. .. 0.1 .. 1.2]
//    |> List.map (fun noiseS ->
//        printfn "noise level: %f" noiseS
//        let idVariance =
//            [1;2;3]
//            |> List.map (fun id ->
//                printfn "insertion: %i" id
//                let averError =
//                    [1 .. 100]
//                    |> List.map (fun x -> 
//                        printfn "iteration %i" x
//                        let noisySignal = synthDeepData noiseS id
//                        let fastTreeNoise = applySSN noisySignal 40 |> Tree.filterLeaves |> Array.rev
////                        let errorPure =
////                            (fastTreeNoise |> List.find (fun xl -> xl |> List.exists (fun p -> p.ID=id)) |> List.length)-1
////                        printfn "purity error: %i" errorPure
//                        let error = 
//                            fastTreeNoise.Length - id - 1//(fastTreeNoise |> List.findIndex (fun xl -> xl |> List.exists (fun p -> p.ID=id))) - 2
//                        printfn "error: %i" error
//                        float (abs(error))  //+ errorPure
//                    )
//                (id, averError) 
//            )   
//        (noiseS, idVariance)
//    )


//let data = (synthDeepData 0.2 3)
//drawKinetik data [|1.; 2.; 3.; 4.; 5.; 6.|] 
//let fastTreeNoise' = applySSN data 10 |> Tree.filterLeaves

//let id1 = checkStat2 |> List.map (fun (n,idList) -> idList.[0] |> snd |> List.map (fun x -> (n-0.02,x/3.))) |> List.concat |> List.unzip
//let id2 = checkStat2 |> List.map (fun (n,idList) -> idList.[1] |> snd |> List.map (fun x -> (n,x/2.))) |> List.concat |> List.unzip
//let id3 = checkStat2 |> List.map (fun (n,idList) -> idList.[2] |> snd |> List.map (fun x -> (n+0.02,x))) |> List.concat |> List.unzip

//let line1 = checkStat2 |> List.map (fun (n,idList) -> (n-0.02, idList.[0] |> snd |> List.averageBy (fun x -> (x/3.)) ))
//let line2 = checkStat2 |> List.map (fun (n,idList) -> (n, idList.[1] |> snd |> List.averageBy (fun x -> (x/2.)) ))
//let line3 = checkStat2 |> List.map (fun (n,idList) -> (n+0.02, idList.[2] |> snd |> List.averageBy (fun x -> (x/1.)) ))

//[Chart.BoxPlot((fst id1),(snd id1),"insertion at level 1");
//Chart.BoxPlot((fst id2),(snd id2),"insertion at level 2");
//Chart.BoxPlot((fst id3),(snd id3),"insertion at level 3");
//Chart.Line(line1,"average, level 1");
//Chart.Line(line2,"average, level 2");
//Chart.Line(line3,"average, level 3");
//]
//|> Chart.Combine
//|> Chart.Show
                        
//// stats 3 - minority revealing

//let dataSynF noise size d = // d - deviation
//    [|
//    [|1 .. (size-d)|] |> Array.map (fun x -> fillProtein [|"1";"1"|] 1. noise);     // major
//    [|1 .. d|] |> Array.map (fun x -> fillProtein [|"1";"1"|] -1. noise);           // minor
//    [|1 .. (size-d)|] |> Array.map (fun x -> fillProtein [|"1";"2"|] -1. noise);    // major
//    [|1 .. d|] |> Array.map (fun x -> fillProtein [|"1";"2"|] 1. noise);            // minor
//    |]
//    |> Array.concat
//    |> Array.mapi (fun i p -> {p with ID=i; ProteinL= Array.init 2 (fun i -> if i=0 then p.ProteinL.[0] else (string i))})

//let dataSyn = dataSynF 0.3 10 4
//let tree = applySSN dataSyn (dataSyn.Length*2)  

//drawKinetikTitle tree ["1";"1";"1"] [|1.; 2.; 3.; 4.; 5.; 6.|]
//sendToGephiFromTreeParam tree

//let checkStat3 =
//    [0. .. 0.1 .. 1.5]
//    |> List.map (fun noiseS ->
//        printfn "noise level: %f" noiseS
//        let idVariance =
//            [1;2;3;4;5]
//            |> List.map (fun d ->
//                printfn "insertion: %i" d
//                let averError =
//                    [1 .. 10]
//                    |> List.map (fun x -> 
//                        printfn "iteration %i" x
//                        let noisySignal = dataSynF noiseS 10 d
//                        let fastTreeNoise = applySSN noisySignal 40 |> Tree.filterLeaves
//                        let minS = 4
//                        let maxS = 20
//                        let error1 = 
//                            fastTreeNoise.Length - minS
//                        let error2 =
//                            fastTreeNoise 
//                            |> Array.map (fun leaf -> 
//                                leaf 
//                                |> Array.groupBy (fun i -> i.ProteinL.[0]) 
//                                |> Array.map (fun (pattern,a) -> (a.Length)) 
//                                |> (fun arr -> 
//                                    if arr.Length=1 then 0 
//                                    else (arr |> Array.min)))
//                            |> Array.sum
//                        let errorRatio = 
//                            (float (error1 + error2))/(float (maxS-minS+error2))
//                        printfn "type I error: %i" error1
//                        printfn "type II error: %i" error2
//                        printfn "error: %f" errorRatio
//                        (float (abs(errorRatio)),error1,error2) 
//                    )
//                (d, averError) 
//            )   
//        (noiseS, idVariance)
//    )

//let d1 = checkStat3 |> List.map (fun (n,idList) -> idList.[0] |> snd |> List.map (fun (aver,er1,er2) -> (n-0.02,aver))) |> List.concat |> List.unzip
//let d2 = checkStat3 |> List.map (fun (n,idList) -> idList.[1] |> snd |> List.map (fun (aver,er1,er2) -> (n-0.01,aver))) |> List.concat |> List.unzip
//let d3 = checkStat3 |> List.map (fun (n,idList) -> idList.[2] |> snd |> List.map (fun (aver,er1,er2) -> (n,aver)))      |> List.concat |> List.unzip
//let d4 = checkStat3 |> List.map (fun (n,idList) -> idList.[3] |> snd |> List.map (fun (aver,er1,er2) -> (n+0.01,aver))) |> List.concat |> List.unzip
//let d5 = checkStat3 |> List.map (fun (n,idList) -> idList.[4] |> snd |> List.map (fun (aver,er1,er2) -> (n+0.02,aver))) |> List.concat |> List.unzip

//let dAverRatio = checkStat3 |> List.map (fun (n,idList) -> idList |> List.map (fun (d,xL) -> xL |> List.map (fun (aver,er1,er2) -> (n,aver))) |> List.concat) |> List.concat |> List.unzip
//let dAverError1 = checkStat3 |> List.map (fun (n,idList) -> idList |> List.map (fun (d,xL) -> xL |> List.map (fun (aver,er1,er2) -> (n,er1))) |> List.concat) |> List.concat |> List.unzip
//let dAverError2 = checkStat3 |> List.map (fun (n,idList) -> idList |> List.map (fun (d,xL) -> xL |> List.map (fun (aver,er1,er2) -> (n,er2))) |> List.concat) |> List.concat |> List.unzip

//let lined1 = checkStat3 |> List.map (fun (n,idList) -> (n-0.02, idList.[0] |> snd |> List.averageBy (fun (aver,er1,er2) -> aver ) ))
//let lined2 = checkStat3 |> List.map (fun (n,idList) -> (n-0.01, idList.[1] |> snd |> List.averageBy (fun (aver,er1,er2) -> aver ) ))
//let lined3 = checkStat3 |> List.map (fun (n,idList) -> (n,      idList.[2] |> snd |> List.averageBy (fun (aver,er1,er2) -> aver ) ))
//let lined4 = checkStat3 |> List.map (fun (n,idList) -> (n+0.01, idList.[3] |> snd |> List.averageBy (fun (aver,er1,er2) -> aver ) ))
//let lined5 = checkStat3 |> List.map (fun (n,idList) -> (n+0.02, idList.[4] |> snd |> List.averageBy (fun (aver,er1,er2) -> aver ) ))

//let linedAverRatio = checkStat3 |> List.map (fun (n,idList) -> (n, idList |> List.map (snd >> List.averageBy (fun (aver,er1,er2) -> aver )) |> List.average))// |> List.averageBy) ))
//let linedAverErrors = checkStat3 |> List.map (fun (n,idList) -> (n, idList |> List.map (snd >> List.averageBy (fun (aver,er1,er2) -> float (er1+er2) )) |> List.average))// |> List.averageBy) ))

//Chart.StackedBar([1,1;1,2;2,2;2,3]) |> Chart.Show  

//[Chart.BoxPlot((fst d1),(snd d1),"ratio 1/10");
//Chart.BoxPlot((fst d2),(snd d2),"ratio 2/10");
//Chart.BoxPlot((fst d3),(snd d3),"ratio 3/10");
//Chart.BoxPlot((fst d4),(snd d4),"ratio 4/10");
//Chart.BoxPlot((fst d5),(snd d5),"ratio 5/10");
//Chart.Line(lined1,"average, 1/10");
//Chart.Line(lined2,"average, 2/10");
//Chart.Line(lined3,"average, 3/10");
//Chart.Line(lined4,"average, 4/10");
//Chart.Line(lined5,"average, 5/10");
//]
//|> Chart.Combine
//|> Chart.Show  

//// or

//[Chart.BoxPlot((fst dAverRatio),(snd dAverRatio),"all ratios together");   
//Chart.Line(linedAverRatio,"average");         
//]
//|> Chart.Combine
//|> Chart.Show   

//// or 

//[Chart.BoxPlot((fst dAverError1),(snd dAverError1),"type I error in all ratios");  
//Chart.BoxPlot((fst dAverError2),(snd dAverError2),"type II error in all ratios");   
//Chart.Line(linedAverErrors,"average");         
//]
//|> Chart.Combine
//|> Chart.Show   

//// draw all ratios separately
//let noise = [0. .. 0.1 .. 1.5]
//let r1 = checkStat3 |> List.map (fun (n,idList) -> idList.[0] |> snd |> List.sumBy (fun (aver,er1,er2) -> er1+er2))
//let r2 = checkStat3 |> List.map (fun (n,idList) -> idList.[1] |> snd |> List.sumBy (fun (aver,er1,er2) -> er1+er2)) 
//let r3 = checkStat3 |> List.map (fun (n,idList) -> idList.[2] |> snd |> List.sumBy (fun (aver,er1,er2) -> er1+er2)) 
//let r4 = checkStat3 |> List.map (fun (n,idList) -> idList.[3] |> snd |> List.sumBy (fun (aver,er1,er2) -> er1+er2)) 
//let r5 = checkStat3 |> List.map (fun (n,idList) -> idList.[4] |> snd |> List.sumBy (fun (aver,er1,er2) -> er1+er2)) 

//let r1a = checkStat3 |> List.map (fun (n,idList) -> idList.[0] |> snd |> List.averageBy (fun (aver,er1,er2) -> aver))
//let r2a = checkStat3 |> List.map (fun (n,idList) -> idList.[1] |> snd |> List.averageBy (fun (aver,er1,er2) -> aver)) 
//let r3a = checkStat3 |> List.map (fun (n,idList) -> idList.[2] |> snd |> List.averageBy (fun (aver,er1,er2) -> aver)) 
//let r4a = checkStat3 |> List.map (fun (n,idList) -> idList.[3] |> snd |> List.averageBy (fun (aver,er1,er2) -> aver)) 
//let r5a = checkStat3 |> List.map (fun (n,idList) -> idList.[4] |> snd |> List.averageBy (fun (aver,er1,er2) -> aver)) 


//[
//Chart.StackedColumn(noise,r1,"ratio=0.1", Color = Plots.colorBlue);
//Chart.StackedColumn(noise,r2,"ratio=0.2", Color = Plots.colorOrange);
//Chart.StackedColumn(noise,r3,"ratio=0.3", Color = Plots.colorGreen);
//Chart.StackedColumn(noise,r4,"ratio=0.4", Color = Plots.colorYellow);
//Chart.StackedColumn(noise,r5,"ratio=0.5", Color = Plots.colorGray);
//]
//|> Chart.Combine
//|> Chart.Show   

//[
//Chart.Column(noise,r1,"ratio=0.1", Color = Plots.colorBlue);
//Chart.Column(noise,r2,"ratio=0.2", Color = Plots.colorOrange);
//Chart.Column(noise,r3,"ratio=0.3", Color = Plots.colorGreen);
//Chart.Column(noise,r4,"ratio=0.4", Color = Plots.colorYellow);
//Chart.Column(noise,r5,"ratio=0.5", Color = Plots.colorGray);
//]
//|> Chart.Stack( 1, 0.1)
//|> Chart.Show

//[
//Chart.Column(noise,r1a,"ratio=0.1", Color = Plots.colorBlue);
//Chart.Column(noise,r2a,"ratio=0.2", Color = Plots.colorOrange);
//Chart.Column(noise,r3a,"ratio=0.3", Color = Plots.colorGreen);
//Chart.Column(noise,r4a,"ratio=0.4", Color = Plots.colorYellow);
//Chart.Column(noise,r5a,"ratio=0.5", Color = Plots.colorGray);
//]
//|> Chart.Stack( 1, 0.1)
//|> Chart.withSize (650., 400.)
//|> Chart.Show

//// Matrix Correlation Assesment:

// complete random dataset of 100 vectors with size=6
let generator n =
    let r = MathNet.Numerics.Distributions.Normal(0., 1.)
    [| for x in [|1 .. 1 .. n|] -> (r.Sample())|] |> General.zScoreTransform

let synData = 
    [|for a in [0 .. 1 .. 29] -> 
        {
        Types.ID=a;
        Types.ProteinL = [|sprintf "%i" a|];
        Types.BinL = [||];
        Types.OriginalBin = [||];
        Types.dataL = generator 3}|]

let matrix = 
    synData
    |> distMatrixWeightedOf distanceMatrixWeighted None

let time = [|1. .. 1. .. 3.|]

let lines (ids: int list) =
    ids
    |> List.map (fun i -> (synData.[i]))
    

let compare2 (dir: int list list) (tes: int list list) tesTitle = 
    let er1 = List.fold2 (fun acc l1 l2 -> acc + abs ((l1 |> List.length) - (l2 |> List.length))) 0 dir tes
    printfn "length diff for %s = %i" tesTitle er1


// find out the best k from direct clustering and apply silhouete criteria
//let clusterSchemes =
//    [|2 .. 1 .. 99|]
//    |> Array.map (fun k -> 
//        synData
//        |> Array.map (fun protein -> protein.dataL)
//        |> ML.Unsupervised.HierarchicalClustering.generate (weightedEuclidean None) (ML.Unsupervised.HierarchicalClustering.Linker.centroidLwLinker)
//        |> ML.Unsupervised.HierarchicalClustering.cutHClust k
//        |> List.map (List.map (fun i -> synData.[ML.Unsupervised.HierarchicalClustering.getClusterId i]))
//    )

let rv2 (x': matrix) (y': matrix) =
    let x =
        x'
        |> Matrix.toJaggedArray 
        |> Array.map (fun ar -> Array.append [|0.000001|] ar) 
        |> MatrixTopLevelOperators.matrix 
        |> Matrix.transpose
    let y =
        y'
        |> Matrix.toJaggedArray 
        |> Array.map (fun ar -> Array.append [|0.000001|] ar) 
        |> MatrixTopLevelOperators.matrix 
        |> Matrix.transpose
    
    let xTilda = x*x.Transpose - (Matrix.diag ((x*x.Transpose).Diagonal))
    let yTilda = y*y.Transpose - (Matrix.diag ((y*y.Transpose).Diagonal))
    let xy = (xTilda.ToVector()).Transpose*(yTilda.ToVector())
    let xx = (xTilda.ToVector()).Transpose*(xTilda.ToVector())
    let yy = (yTilda.ToVector()).Transpose*(yTilda.ToVector())
    xy/sqrt(xx*yy)

let pairwiseCorrMax (x:matrix) (y:matrix) =
    let xN = x.Dimensions |> fst
    let yN = y.Dimensions |> fst
    let m = Array2D.create xN yN 0.
    for rowI in 0..(xN-1) do
        for colI in 0..(yN-1) do
            let tmp = weightedEuclidean None (x.Row rowI) (y.Row colI) 
            m.[rowI,colI] <- tmp
    m 
    |> Array2D.array2D_to_seq 
    |> Seq.max

let pairwiseCorrSum (x:matrix) (y:matrix) =
    let xN = x.Dimensions |> fst
    let yN = y.Dimensions |> fst
    let m = Array2D.create xN yN 0.
    for rowI in 0..(xN-1) do
        for colI in 0..(yN-1) do
            let tmp = weightedEuclidean None (x.Row rowI) (y.Row colI) 
            m.[rowI,colI] <- tmp
    m 
    |> Array2D.array2D_to_seq 
    |> Seq.sum

let pairwiseCorrAverage (x:matrix) (y:matrix) =
    let xN = x.Dimensions |> fst
    let yN = y.Dimensions |> fst
    let m = Array2D.create xN yN 0.
    for rowI in 0..(xN-1) do
        for colI in 0..(yN-1) do
            let tmp = weightedEuclidean None (x.Row rowI) (y.Row colI) 
            m.[rowI,colI] <- tmp
    m 
    |> Array2D.array2D_to_seq 
    |> Seq.average

let svdDCorr1 (x:matrix) (y:matrix) =
    let svdF mat =
        mat |> Matrix.toJaggedArray |> DenseMatrix.OfRowArrays |> fun m ->
            let svd = m.Svd()
            let s = svd.S.[0]
            (svd.VT.Row 0).Map (fun i -> i*s) 
    let xS = svdF x
    let yS = svdF y
    euclidean xS yS

let svdDCorr2 (x:matrix) (y:matrix) =
    let svd2F mat =
        mat 
        |> Matrix.toJaggedArray 
        |> DenseMatrix.OfRowArrays
        |> (fun m ->
            let svd = m.Svd()
            let rowN = if m.RowCount < 2 then 1 else 2
            svd.VT.ToRowArrays()
            |> Array.truncate rowN
            |> Array.mapi (fun id v ->
                let s = svd.S.[id]
                v |> Array.map (fun i -> i*s) )
            |> Array.fold (fun acc vx -> vx |> Array.mapi (fun i vv -> vv + acc.[i])) (Array.zeroCreate (m.ColumnCount)) 
                     )
    let xS = svd2F x
    let yS = svd2F y
    euclidean xS yS

let svdDCorrAll (x:matrix) (y:matrix) =
    let svdAllF mat =
        mat 
        |> Matrix.toJaggedArray 
        |> DenseMatrix.OfRowArrays 
        |> (fun m ->
            let svd = m.Svd()
            svd.VT.ToRowArrays()
            |> Array.truncate m.RowCount
            |> Array.mapi (fun id v ->
                let s = svd.S.[id]
                v |> Array.map (fun i -> i*s) )
            |> Array.fold (fun acc vx -> vx |> Array.mapi (fun i vv -> vv + acc.[i])) (Array.zeroCreate (m.ColumnCount)) 
                     )
                                                    
    let xS = svdAllF x
    let yS = svdAllF y
    euclidean xS yS

let svdRCorr1 (x:matrix) (y:matrix) =
    let svdF mat =
        mat |> Matrix.toJaggedArray |> DenseMatrix.OfRowArrays |> fun m ->
            let svd = m.Svd()
            let s = svd.S.[0]
            (svd.VT.Row 0).Map (fun i -> i*s) 
    let xS = svdF x
    let yS = svdF y
    pearson xS yS

let svdRCorr2 (x:matrix) (y:matrix) =
    let svd2F mat =
        mat 
        |> Matrix.toJaggedArray 
        |> DenseMatrix.OfRowArrays 
        |> (fun m ->
            let svd = m.Svd()
            let rowN = if m.RowCount < 2 then 1 else 2
            svd.VT.ToRowArrays()
            |> Array.truncate rowN
            |> Array.mapi (fun id v ->
                let s = svd.S.[id]
                v |> Array.map (fun i -> i*s) )
            |> Array.fold (fun acc vx -> vx |> Array.mapi (fun i vv -> vv + acc.[i])) (Array.zeroCreate (m.ColumnCount)) 
                     )
    let xS = svd2F x
    let yS = svd2F y
    pearson xS yS

let svdRCorrAll (x:matrix) (y:matrix) =
    let svdAllF mat =
        mat 
        |> Matrix.toJaggedArray 
        |> DenseMatrix.OfRowArrays 
        |> (fun m ->
            let svd = m.Svd()
            svd.VT.ToRowArrays()
            |> Array.truncate m.RowCount
            |> Array.mapi (fun id v ->
                let s = svd.S.[id]
                v |> Array.map (fun i -> i*s) )
            |> Array.fold (fun acc vx -> vx |> Array.mapi (fun i vv -> vv + acc.[i])) (Array.zeroCreate (m.ColumnCount)) 
                     )
                                                    
    let xS = svdAllF x
    let yS = svdAllF y
    pearson xS yS

// set k_pre > k_end
let k_pre = 12
let k_end = 4

// direct clustering with k=k_end
let end_result =
    synData
    |> Array.map (fun protein -> protein.dataL)
    |> ML.Unsupervised.HierarchicalClustering.generate (weightedEuclidean None) (ML.Unsupervised.HierarchicalClustering.Linker.completeLwLinker)
    |> ML.Unsupervised.HierarchicalClustering.cutHClust k_end
    |> List.map (List.map (fun i -> synData.[ML.Unsupervised.HierarchicalClustering.getClusterId i]))
    |> List.map (fun cluster -> 
        let bin = cluster |> List.map (fun p -> p.ID) |> List.sort
        let matrixC = MatrixTopLevelOperators.matrix (cluster |> List.map (fun p -> p.dataL) |> List.toArray)
        (bin,matrixC))

// make pre-clustering with k=k_pre
let pre_clusters = 
    synData
    |> Array.map (fun protein -> protein.dataL)
    |> ML.Unsupervised.HierarchicalClustering.generate (weightedEuclidean None) (ML.Unsupervised.HierarchicalClustering.Linker.completeLwLinker)
    |> ML.Unsupervised.HierarchicalClustering.cutHClust k_pre
    |> List.map (List.map (fun i -> synData.[ML.Unsupervised.HierarchicalClustering.getClusterId i]))
    |> List.map (fun cluster -> 
        let bin = cluster |> List.map (fun p -> p.ID) |> List.sort
        let matrixC = MatrixTopLevelOperators.matrix (cluster |> List.map (fun p -> p.dataL) |> List.toArray)
        (bin,matrixC))

// cluster resulting matrices with matrix correlation as a distance
let matrix_clustering_result fn =
    let clusters nClusters =
        let rvToMap (mapA: ((int list) * matrix))  (mapB: ((int list) * matrix)) = 
            let mA = 
                mapA 
                |> snd 
            let mB = 
                mapB 
                |> snd 
            fn mA mB 
// rv2 // pairwiseCorrMax // pairwiseCorrSum // pairwiseCorrAverage // svdDCorr1 // svdDCorr2 // svdDCorrAll // svdRCorr1 // svdRCorr2 // svdRCorrAll
        pre_clusters
        |> ML.Unsupervised.HierarchicalClustering.generate rvToMap (ML.Unsupervised.HierarchicalClustering.Linker.completeLwLinker) // here!
        |> ML.Unsupervised.HierarchicalClustering.cutHClust nClusters
        |> List.map (List.map (fun i -> pre_clusters.[ML.Unsupervised.HierarchicalClustering.getClusterId i]))
    clusters k_end
    |> List.map (fun list ->
        let binName = list |> List.map (fun (bin,p) -> bin) |> List.concat |> List.sort
        let items = list |> List.map snd
        (binName,items))

// compare two results (clusters content)

let direct = end_result |> List.map fst |> List.sortBy (fun x -> x.Length)

let twostepRV2 = matrix_clustering_result rv2 |> List.map fst |> List.sortBy (fun x -> x.Length)
let twostepPCMax = matrix_clustering_result pairwiseCorrMax |> List.map fst |> List.sortBy (fun x -> x.Length)
let twostepPCSum = matrix_clustering_result pairwiseCorrSum |> List.map fst |> List.sortBy (fun x -> x.Length)
let twostepPCAve = matrix_clustering_result pairwiseCorrAverage |> List.map fst |> List.sortBy (fun x -> x.Length)
let twostepSVDD1 = matrix_clustering_result svdDCorr1 |> List.map fst |> List.sortBy (fun x -> x.Length)
let twostepSVDD2 = matrix_clustering_result svdDCorr2 |> List.map fst |> List.sortBy (fun x -> x.Length)
let twostepSVDDA = matrix_clustering_result svdDCorrAll |> List.map fst |> List.sortBy (fun x -> x.Length)
let twostepSVDR1 = matrix_clustering_result svdRCorr1 |> List.map fst |> List.sortBy (fun x -> x.Length)
let twostepSVDR2 = matrix_clustering_result svdRCorr2 |> List.map fst |> List.sortBy (fun x -> x.Length)
let twostepSVDRA = matrix_clustering_result svdRCorrAll |> List.map fst |> List.sortBy (fun x -> x.Length)

compare2 direct twostepRV2 "twostepRV2"
compare2 direct twostepPCMax "twostepPCMax"
compare2 direct twostepPCSum "twostepPCSum"
compare2 direct twostepPCAve "twostepPCAve"
compare2 direct twostepSVDD1 "twostepSVDD1"
compare2 direct twostepSVDD2 "twostepSVDD2"
compare2 direct twostepSVDDA "twostepSVDDA"
compare2 direct twostepSVDR1 "twostepSVDR1"
compare2 direct twostepSVDR2 "twostepSVDR2"
compare2 direct twostepSVDRA "twostepSVDRA"

// don't coincide!!!!!!!!!

direct
|> List.mapi (fun c list ->
    let col =
        match c with
        |0 -> colorBlue
        |1 -> colorGray
        |2 -> colorGreen
        |_ -> colorOrange
    list |> lines |> List.map (fun d -> Chart.Line(time, d.dataL, sprintf "%i" d.ID, Color=col) ) |> Chart.Combine
    )
|> Chart.Combine
|> Chart.withTitle "Single Linker: direct"
|> Chart.Show

//
//// K-Mean testing setup
//

let centroidRandom (input: float [] array) (k: int) : (float [] []) =
    let r = new System.Random() 
    IterativeClustering.randomCentroids r input k

// direct clustering with k=k_end
let end_result_kmean =
    let data = 
        synData
        |> Array.map (fun protein -> protein.dataL)
    let clusters = 
        let c1 = ML.Unsupervised.IterativeClustering.kmeans euclidean (centroidRandom) data k_end
        let x1 = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c1
        [|1 .. 20|]
        |> Array.fold (fun (disp,best) x -> 
            let c = ML.Unsupervised.IterativeClustering.kmeans euclidean (centroidRandom) data k_end
            let x = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c
            if x<disp then
                (x,c)
            else
                (disp,best) ) (x1,c1)
        |> snd

    synData
    |> Array.map (fun item -> (clusters.Classifier item.dataL |> fst), item)
    |> Array.groupBy fst
    |> Array.map (fun (c,group) -> 
        group 
        |> Array.map snd
        |> (fun cluster -> 
                            let bin = cluster |> Array.map (fun p -> p.ID) |> Array.sort |> List.ofArray
                            let matrixC = MatrixTopLevelOperators.matrix (cluster |> Array.map (fun p -> p.dataL) )
                            (bin,matrixC)))
    |> List.ofArray 
 
let direct_kmean = end_result_kmean |> List.map fst |> List.sortBy (fun x -> x.Length)


// make pre-clustering with k=k_pre
let pre_clusters_kmean = 
    let data = 
        synData
        |> Array.map (fun protein -> protein.dataL)
    let clusters = 
        let c1 = ML.Unsupervised.IterativeClustering.kmeans euclidean (ML.Unsupervised.IterativeClustering.intitCVMAX) data k_pre
        let x1 = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c1
        [|1 .. 20|]
        |> Array.fold (fun (disp,best) x -> 
            let c = ML.Unsupervised.IterativeClustering.kmeans euclidean (ML.Unsupervised.IterativeClustering.intitCVMAX) data k_pre
            let x = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c
            if x<disp then
                (x,c)
            else
                (disp,best) ) (x1,c1)
        |> snd
    synData
    |> Array.map (fun item -> (clusters.Classifier item.dataL |> fst), item)
    |> Array.groupBy fst
    |> Array.map (fun (cluster,group) -> 
        group 
        |> Array.map snd
        |> fun cluster -> 
            let bin = cluster |> Array.map (fun p -> p.ID) |> Array.sort |> List.ofArray
            let matrixC = MatrixTopLevelOperators.matrix (cluster |> Array.map (fun p -> p.dataL) )
            (bin,matrixC))
    |> List.ofArray

//// Pairwise

let distMatrix (matrixA: (int list) * matrix) (matrixB: (int list) * matrix) =
    let mA = snd matrixA
    let mB = snd matrixB
    pairwiseCorrAverage mA mB

let centroidFactory (input: ((int list)*matrix) array) (k: int) : (((int list)*matrix) []) =
    let r = new System.Random() 
    IterativeClustering.randomCentroids r input k


let intiCgroups (input: ((int list)*matrix) array) k : (((int list)*matrix) []) =
        let dmatrix = input |> Array.map (fun (bin,data) -> (data |> Matrix.meanColumnWise)) |> MatrixTopLevelOperators.matrix
        let cvmax = // find a feature with the biggest variance and return the (row Number, values of the feature), sorted by the values 
            dmatrix
            |> Matrix.Generic.enumerateColumnWise Seq.var
            |> Seq.zip (Matrix.Generic.enumerateColumnWise id dmatrix)
            |> Seq.maxBy snd
            |> fst
            |> Seq.mapi (fun rowI value -> (rowI,value)) 
            |> Seq.toArray 
            |> Array.sortBy snd
                    
        if cvmax.Length < k then failwithf "Number of data points must be at least %i" k        
//        //
//        let intDivide a b = 
//            int (System.Math.Floor((float a) / (float b)))
    
        let chunkSize = cvmax.Length / k
        let midChunk  = chunkSize / 2
        [ for i=1 to k do
            let index = 
                match (chunkSize * i) with
                | x when x < cvmax.Length -> x - midChunk
                | x                       -> chunkSize * (i - 1) + ((cvmax.Length - chunkSize * (i - 1)) / 2)
            //printfn "Array.lenght = %i and index = %i" cvmax.Length (index-1)
            yield cvmax.[index-1] |> fst]
        |> Seq.map (fun rowI -> input.[rowI])
        |> Seq.toArray

// Recompute Centroid as average of given sample (for kmeans)
let updateCentroid (current: (int list) * matrix) (sample: ((int list) * matrix) []) = // rewrite it in matrix!
    let size = sample.Length
    match size with
    | 0 -> current
    | _ ->
        ([], 
            sample
            |> Array.map (fun (_,x) -> x.ToArray2D() |> Array2D.toJaggedArray)
            |> Array.concat
            |> fun p -> MatrixTopLevelOperators.matrix p)

let matrix_kmean_pd =    
    let clusters = 
        let c1 = ML.Unsupervised.IterativeClustering.compute distMatrix (intiCgroups) updateCentroid (pre_clusters_kmean |> Array.ofList) k_end
        let x1 = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c1
        [|1 .. 20|]
        |> Array.fold (fun (disp,best) x -> 
            let c = ML.Unsupervised.IterativeClustering.compute distMatrix (intiCgroups) updateCentroid (pre_clusters_kmean |> Array.ofList) k_end
            let x = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c
            if x<disp then
                (x,c)
            else
                (disp,best) ) (x1,c1)
        |> snd

    pre_clusters_kmean
    |> List.map (fun list -> (clusters.Classifier list |> fst),list )
    |> List.groupBy (fst)
    |> List.map (fun (cID,list) ->
        let binName = list |> List.map (fun (cID,(bin,p)) -> bin) |> List.concat |> List.sort
        let items = list |> List.map (fun (cID,(bin,p)) -> p)
        (binName,items))
    |> List.map fst |> List.sortBy (fun x -> x.Length)

compare2 direct matrix_kmean_pd "matrix_kmean_pd"

matrix_kmean_pd
|> List.mapi (fun c list ->
    let col =
        match c with
        |0 -> colorBlue
        |1 -> colorGray
        |2 -> colorGreen
        |_ -> colorOrange
    list |> lines |> List.map (fun d -> Chart.Line(time, d.dataL, sprintf "%i" d.ID, Color=col) ) |> Chart.Combine
    )
|> Chart.Combine
|> Chart.withTitle "K-Mean: matrix_kmean_pd"
|> Chart.Show


//// SVD-based

// Given a distance, centroid factory and
// centroid aggregation function, identify
// the k centroids of a dataset
let computeExp (dist: Distance<'a>) 
            (factory: FSharp.Stats.ML.Unsupervised.IterativeClustering.CentroidsFactory<'a>) 
            (aggregator: FSharp.Stats.ML.Unsupervised.IterativeClustering.ToCentroid<'a>)
            (dataset: 'a array) 
            k =
    
    let closest (dist: Distance<'a>) centroids (obs: 'a) =
        centroids
        |> Array.mapi (fun i c -> (i, dist c obs)) 
        |> Array.minBy (fun (i, d) -> d)

    // Recursively update Centroids and
    // the assignment of observations to Centroids
    let rec update (centroids, assignment) =
        // Assign each point to the closest centroid
        let next = 
            dataset 
            |> Array.map (fun obs -> closest dist centroids obs)
            //|> Seq.toList
        // Check if any assignment changed
        let change =
            match assignment with
            | Some(previous) -> 
                Array.zip previous next    
                |> Array.exists (fun ((i, _), (j, _)) -> not (i = j))
            | None -> true // initially we have no assignment
        printfn "change is there: %b" change
        if change 
        then 
            // Update each Centroid position:
            // extract cluster of points assigned to each Centroid
            // and compute the new Centroid by aggregating cluster
            let updatedCentroids =
                let assignedDataset = Array.zip dataset next
                centroids 
                |> Array.mapi (fun i centroid -> 
                    assignedDataset 
                    |> Array.filter (fun (_, (ci, _)) -> ci = i)
                    |> Array.map (fun (obs, _) -> obs)
                    |> aggregator centroid)
            // Perform another round of updates
            update (updatedCentroids, Some(next))
        // No assignment changed, we are done
        else (centroids, next)

    let initialCentroids = factory dataset k
//        let centroids = update (initialCentroids, None) |> fst |> Seq.toList        
//        let classifier = fun datapoint -> 
//            centroids 
//            |> List.minBy (fun centroid -> dist centroid datapoint)        
    let centroids,closestDistances = update (initialCentroids, None)        
    let lCentroids = Seq.zip [1..k] centroids |> Seq.toArray
    let classifier = fun datapoint -> 
        lCentroids 
        |> Array.minBy (fun centroid -> dist (snd centroid) datapoint)
    FSharp.Stats.ML.Unsupervised.IterativeClustering.createKClusteringResult lCentroids classifier closestDistances dist        


let svdFm mat =
    mat 
    |> Matrix.toJaggedArray 
    |> DenseMatrix.OfRowArrays 
    |> (fun m ->
        let svd = m.Svd()
        let sN = if m.RowCount < 2 then 1 else 2
        svd.VT.ToRowArrays()
        |> Array.truncate sN
        |> Array.mapi (fun id v ->
            let s = svd.S.[id]
            v |> Array.map (fun i -> i*s) ))
    |> MatrixTopLevelOperators.matrix

let svdF mat =
    mat 
    |> Matrix.toJaggedArray 
    |> DenseMatrix.OfRowArrays 
    |> (fun m ->
        let svd = m.Svd()
        let sN = if m.RowCount < 2 then 1 else 2
        svd.VT.ToRowArrays()
        |> Array.truncate sN
        |> Array.mapi (fun id v ->
            let s = svd.S.[id]
            v |> Array.map (fun i -> i*s) )
        |> Array.fold (fun acc vx -> vx |> Array.mapi (fun i vv -> vv + acc.[i])) (Array.zeroCreate (m.ColumnCount)) 
        )

let distMatrixSVD (matrixA: (int list) * matrix) (matrixB: (int list) * matrix) =
    let mA = snd matrixA |> svdF
    let mB = snd matrixB |> svdF
    euclidean mA mB // pearson

let centroidFactorySVD (input: ((int list)*matrix) array) (k: int) : (((int list)*matrix) []) =
    let r = new System.Random() 
    IterativeClustering.randomCentroids r input k
    |> Array.map (fun (b,x) -> (b, x |> svdFm))

let updateCentroidSVD (current: (int list) * matrix) (sample: ((int list) * matrix) []) = 
    let size = sample.Length
    match size with
    | 0 -> current
    | _ ->
        ([], sample
        |> Array.map (fun (_,x) -> x.ToArray2D() |> Array2D.toJaggedArray)
        |> Array.concat
        |> MatrixTopLevelOperators.matrix
        |> svdFm 
        )

let matrix_kmean_svd =
    let clusters = 
        let c1 = computeExp distMatrixSVD (centroidFactorySVD) updateCentroidSVD (pre_clusters_kmean |> Array.ofList) k_end
        let x1 = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c1
        [|1 .. 20|]
        |> Array.fold (fun (disp,best) x -> 
            let c = computeExp distMatrixSVD (centroidFactorySVD) updateCentroidSVD (pre_clusters_kmean |> Array.ofList) k_end
            let x = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c
            if x<disp then
                (x,c)
            else
                (disp,best) ) (x1,c1)
        |> snd

    pre_clusters_kmean
    |> List.map (fun list -> (clusters.Classifier list |> fst),list )
    |> List.groupBy (fst)
    |> List.map (fun (cID,list) ->
        let binName = list |> List.map (fun (cID,(bin,p)) -> bin) |> List.concat |> List.sort
        let items = list |> List.map (fun (cID,(bin,p)) -> p)
        (binName,items))
    |> List.map fst |> List.sortBy (fun x -> x.Length)

compare2 direct matrix_kmean_svd "matrix_kmean_svd"

matrix_kmean_svd
|> List.mapi (fun c list ->
    let col =
        match c with
        |0 -> colorBlue
        |1 -> colorGray
        |2 -> colorGreen
        |_ -> colorOrange
    list |> lines |> List.map (fun d -> Chart.Line(time, d.dataL, sprintf "%i" d.ID, Color=col) ) |> Chart.Combine
    )
|> Chart.Combine
|> Chart.withTitle "K-Mean: SVD-based"
|> Chart.Show


///// K-Mean-Swap testing setup

let sample = synData |> Array.map (fun p -> p.dataL)

let intitCVMAX (sample: float[] array) k =
        let dmatrix = MatrixTopLevelOperators.matrix sample
        let cvmax = // find a feature with the biggest variance and return the (row Number, values of the feature), sorted by the values 
            dmatrix
            |> Matrix.Generic.enumerateColumnWise Seq.var
            |> Seq.zip (Matrix.Generic.enumerateColumnWise id dmatrix)
            |> Seq.maxBy snd
            |> fst
            |> Seq.mapi (fun rowI value -> (rowI,value)) 
            |> Seq.toArray 
            |> Array.sortBy snd
                    
        if cvmax.Length < k then failwithf "Number of data points must be at least %i" k        
//        //
//        let intDivide a b = 
//            int (System.Math.Floor((float a) / (float b)))
    
        let chunkSize = cvmax.Length / k
        let midChunk  = chunkSize / 2
        [ for i=1 to k do
            let index = 
                match (chunkSize * i) with
                | x when x < cvmax.Length -> x - midChunk
                | x                       -> chunkSize * (i - 1) + ((cvmax.Length - chunkSize * (i - 1)) / 2)
            //printfn "Array.lenght = %i and index = %i" cvmax.Length (index-1)
            yield cvmax.[index-1] |> fst]
        |> Seq.map (fun rowI -> dmatrix.Row(rowI).ToArray())
        |> Seq.toArray

let setN = 60
let nRoot = 30
//let depth = 0

let gainFn = (SSN.getStepGainNodeSetnR setN)

let result_kmswap depth power k data =
    
    let parentItems =
        data
        |> Map.toArray
        |> Array.map snd
        |> Array.concat

    let loop (nodeMembers: Types.Item array) dPredSum =
        let dCurrSum = dSumFn (groupIDFn nodeMembers) (groupIDFn nodeMembers) matrix

        /// to calc step from parent to current node
        let stepGain = (gainFn dCurrSum dPredSum nodeMembers.Length nRoot)

        let children = Map.empty
        
        let confGain = confGainFn children
        {
        Member = nodeMembers;
        Children = children
        StepGain = stepGain; 
        ConfGain = (confGain, children |> Map.toList |> List.map fst);
        GroupGain = max stepGain confGain;
        }

    data
    |> FunctionsExp.KMeanSwapFunctions.kmeanSwapShuffleOldsetK power gainFn matrix depth k
    |> Seq.fold (fun (singles,best) i -> 
                        let newNodes = 
                            if singles=Map.empty then
                                i
                                |> Map.fold (fun state key nodes ->
                                                    let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn parentItems) matrix
                                                    state 
                                                    |> Map.add key (loop nodes dPredSum')) (Map.empty)
                            else
                                i
                                |> Map.fold (fun state key nodes ->
                                                    match (singles.TryFind key) with
                                                    | None ->
                                                        let dPredSum' = dSumFn (groupIDFn nodes) (groupIDFn parentItems) matrix
                                                        state 
                                                        |> Map.add key (loop nodes dPredSum')
                                                    | Some x ->
                                                        state 
                                                        |> Map.add key x                                                                   
                                                ) (Map.empty)
                                                                    
                        let best' =
                            if (confGainFn newNodes) > (confGainFn best) then  // compare configuration gains to get the best
                                newNodes  
                            else 
                                best
                        if (singles = Map.empty) then
                            (newNodes, best')
                        else
                            (singles, best')
                                            
                    ) (Map.empty,Map.empty) // here as state should be this singles (first) saved and optimal conf
    |> snd
    |> Map.toList
    |> List.map (fun (s,x) -> 
        (x, s 
            |> String.split '-' |> (fun ar -> if ar.Length>1 then ar.[1 ..] else ar) |> Array.map int |> Array.sort |> Array.toList))
 
//let direct_ideal =
//    let data = synData |> Array.truncate 10 |> Array.map (fun p -> {p with BinL = [|"1"|]; OriginalBin = [|"1"|]})
//    applySSNcombi data setN

//let direct_kmswap_ideal = 
//    direct_ideal.Children
//    |> Map.toList
//    |> List.map (fun x -> 
//          x 
//          |> fst 
//          |> String.replace "p" "" 
//          |> String.split '-' 
//          |> (fun ar -> if ar.Length>1 then ar.[1 ..] else ar) 
//          |> Array.map int 
//          |> Array.sort 
//          |> Array.toList)
//    |> List.sortBy (fun x -> x.Length)

let direct_kmswap = 
    synData 
    |> Array.map (fun p -> (string p.ID, [|{p with BinL = [|"1"|]; OriginalBin = [|"1"|]}|])) 
    |> Map.ofArray 
    |> result_kmswap 0 100 k_end 
    //|> List.sortBy (fun x -> x.Length)

let pre_clusters_kmswap = pre_clusters_kmean

let matrix_kmeanswap_result =
    pre_clusters_kmswap 
    |> List.map (fun (bin,_) -> 
                            let label = String.Join("-",bin)
                            let items = bin |> List.map (fun i -> {synData.[i] with BinL = [|"1";label|]; OriginalBin = [|"1"|]}) |> List.toArray
                            (label,items))
    |> Map.ofList 
    |> result_kmswap 1 50 k_end 
    //|> List.sortBy (fun x -> x.Length)



matrix_kmeanswap_result
|> List.map snd
|> List.sortBy (fun x -> x.Length)
|> List.mapi (fun c list ->
    let col =
        match c with
        |0 -> colorBlue
        |1 -> colorGray
        |2 -> colorGreen
        |3 -> colorOrange
        |4 -> colorYellow
        |5 -> "rgba(230,230,230,1)"
        |6 -> "rgba(20,180,180,1)"
        |7 -> "rgba(180,20,180,1)"
        |8 -> "rgba(180,180,20,1)"
        |9 -> "rgba(20,20,180,1)"
        |10 -> "rgba(180,20,20,1)"
        |_ -> colorOrange
    list |> lines |> List.map (fun d -> Chart.Line(time, d.dataL, sprintf "%i" d.ID, Color=col) ) |> Chart.Combine
    )
|> Chart.Combine
|> Chart.withTitle "K-Mean Swap: matrix_kmeanswap_result"
|> Chart.Show

