//#r @"c:\Users\mikha\source\repos\mathnet-numerics\src\Numerics\bin\Debug\netstandard2.0\MathNet.Numerics.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\netstandard.dll"

//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\BioFSharp.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\BioFSharp.IO.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharpAux.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharpAux.IO.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\Newtonsoft.Json.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharp.Stats.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharp.Plotly.dll"
//#r @"c:\Users\mikha\source\repos\SSN\SSN\bin\Debug\FSharpGephiStreamer.dll"

//#load "Types.fs"
//#load "Functions.fs"
//#load "FunctionsExp.fs"
//#load "GePhi.fs"
//#load "TestData.fs"
//#load "SOM.fs"
//#load "Auxilliary.fs"
//#load "Plots.fs"

//open System 
//open FSharpAux
//open FSharp.Plotly
//open FSharp.Stats

//open Functions
//open Functions.SSN
//open Functions.KMeanSwapFunctions
//open TestData
//open GePhi
//open Types
//open Auxilliary
//open Plots
//open FunctionsExp

//let generateSyn pattern noiseSigma n =
//    let r = MathNet.Numerics.Distributions.Normal(0., noiseSigma)
//    [|1 .. n|] 
//    |> Array.mapFold (fun prev i -> (pattern*log(prev*(float i))),(prev+1.)) 1.
//    |> fst 
//    |> General.zScoreTransform
//    |> Array.map (fun x -> x + (r.Sample()))
//    |> General.zScoreTransform

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

//let pureSignal = 
//    Array.append (synData1Level -1. 0. 2) (synData1Level 1. 0. 2) 
//    |> Array.mapi (fun i p -> {p with ID=i})  

//let noisySignal noise = 
//    Array.append (synData1Level 1. noise 5) (synData1Level -1. noise 5) 
//    |> Array.mapi (fun i p -> {p with ID=i;ProteinL=[|sprintf "%s-%i" p.ProteinL.[0] i|]})

//let data1 = // 11111-22222
//    [|{ID = 0;
//        ProteinL = [|"name-1-0"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.594924196; -0.6799640241; -0.0256084639; 0.4184221454; 0.7876903374; 1.094384202|];};
//    {ID = 1;
//        ProteinL = [|"name-1-1"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.623806497; -0.6918776997; 0.01799425629; 0.485284777; 0.8499604241; 0.9624447397|];};
//    {ID = 2;
//        ProteinL = [|"name-1-2"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.660086232; -0.6015048645; -0.01292323899; 0.5076795715; 0.7051654513; 1.061669313|];};
//    {ID = 3;
//        ProteinL = [|"name-1-3"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.679000203; -0.5919325127; 0.01478577775; 0.4538017553; 0.9110590795; 0.8912861036|];};
//    {ID = 4;
//        ProteinL = [|"name-1-4"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.601316534; -0.737273413; 0.09090209465; 0.4741557857; 0.6788978321; 1.094634234|];};
//    {ID = 5;
//        ProteinL = [|"name-2-5"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.620908382; 0.6052804092; 0.1000292681; -0.4897408503; -0.7310015607; -1.105475648|];};
//    {ID = 6;
//        ProteinL = [|"name-2-6"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.732824906; 0.4496575389; -0.006152173633; -0.358682348; -0.8236386405; -0.9940092832|];};
//    {ID = 7;
//        ProteinL = [|"name-2-7"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.661185627; 0.5405163318; 0.08330491781; -0.4533593937; -0.7449967746; -1.086650709|];};
//    {ID = 8;
//        ProteinL = [|"name-2-8"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.571046142; 0.6659879663; 0.1188755915; -0.4858593321; -0.7234369518; -1.146613416|];};
//    {ID = 9;
//        ProteinL = [|"name-2-9"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.67751769; 0.5858771221; -0.1169415339; -0.4590948617; -0.5317098437; -1.155648573|];}|]

//let data2 = // 11111-22222
//    [|{ID = 0;
//        ProteinL = [|"name-1-0"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.249006804; -1.252379798; 0.307344253; 0.9641723123; 0.4014543971; 0.8284156403|];};
//    {ID = 1;
//        ProteinL = [|"name-1-1"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.904745432; -0.03545858216; 0.235293155; 0.3717654776; 1.04635999; 0.2867853912|];};
//    {ID = 2;
//        ProteinL = [|"name-1-2"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.40088378; -0.4039629479; -0.4222492444; 0.1571901489; 0.5206182222; 1.549287601|];};
//    {ID = 3;
//        ProteinL = [|"name-1-3"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.64537762; -0.3346709741; -0.01776556405; 0.7748367887; 1.256548116; -0.03357074684|];};
//    {ID = 4;
//        ProteinL = [|"name-1-4"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.711570784; -0.1134732279; -0.3248863839; 0.7981286614; 0.2279026571; 1.123899077|];};
//    {ID = 5;
//        ProteinL = [|"name-2-5"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.312453878; 0.6927622691; -1.024884286; -0.5294283502; 0.60074619; -1.051649701|];};
//    {ID = 6;
//        ProteinL = [|"name-2-6"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.041926113; 0.5805275977; 0.6952254848; 0.07572638555; -0.8619569292; -1.531448652|];};
//    {ID = 7;
//        ProteinL = [|"name-2-7"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|0.9182899219; 1.299063406; -0.3361523136; -0.08224203814; -1.504164967; -0.294794009|];};
//    {ID = 8;
//        ProteinL = [|"name-2-8"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.797260723; 0.4243921594; -0.2510013552; -0.2785737992; -0.7522928488; -0.9397848787|];};
//    {ID = 9;
//        ProteinL = [|"name-2-9"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.45649186; 1.009780815; -0.1196038671; -0.7154160467; -0.8499496687; -0.7813030926|];}|]

//let data3 = // 00-11-222222
//    [|{ID = 0;
//        ProteinL = [|"name-0-0"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-0.0636620859; -0.06119196923; -0.01466371581; 0.08523164606; 0.08380272863; 0.07048339629|];};
//    {ID = 1;
//        ProteinL = [|"name-0-1"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-0.0355264886; -0.08083434143; 0.0183349666; -0.03903575135; 0.05908404028; 0.07797757442|];};
//    {ID = 2;
//        ProteinL = [|"name-1-2"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.178029051; -1.187539947; -0.0563792052; 0.4131349135; 0.9320466381; 1.076766652|];};
//    {ID = 3;
//        ProteinL = [|"name-1-3"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.859759663; 0.3637071728; 0.3278055376; -0.3354497542; 0.5808568543; 0.9228398528|];};
//    {ID = 4;
//        ProteinL = [|"name-2-4"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.496163239; 0.6758421078; -0.1958494622; -0.5284020891; -0.03861355314; -1.409140242|];};
//    {ID = 5;
//        ProteinL = [|"name-2-5"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.549159199; 0.4256318769; 0.3750266048; -0.5547714976; -0.4735941401; -1.321452043|];};
//    {ID = 6;
//        ProteinL = [|"name-2-6"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.193397759; 0.9218648743; 0.4772811654; -0.4754283842; -1.183154346; -0.9339610681|];};
//    {ID = 7;
//        ProteinL = [|"name-2-7"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.667690621; 0.4968942602; -0.3091856537; -0.001867291717; -0.6452393139; -1.208292622|];};
//    {ID = 8;
//        ProteinL = [|"name-2-8"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.470419431; 0.4840225401; 0.457031256; -0.4363468966; -1.343758816; -0.6313675145|];};
//    {ID = 9;
//        ProteinL = [|"name-2-9"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.808877991; 0.4811514528; -0.3870354405; -0.6965223446; -0.3446561006; -0.8618155582|];}|]

//let data4 = //// it has to be adjusted
//    [|{ID = 0;
//        ProteinL = [|"name-1-0"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.695517053; 0.1977952906; -0.27546189; 0.5012073488; -0.05323207048; 1.325208374|];};
//    {ID = 1;
//        ProteinL = [|"name-1-1"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.159536014; -1.13053636; 0.1159842631; 0.6251900104; 0.1525542806; 1.39634382|];};
//    {ID = 2;
//        ProteinL = [|"name-1-2"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.518924248; -0.9475592846; 0.2818354427; 0.4399212567; 0.871192729; 0.8735341041|];};
//    {ID = 3;
//        ProteinL = [|"name-1-3"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.918262683; 0.09206142105; 0.1295746959; 0.1813512622; 0.5188014429; 0.9964738607|];};
//    {ID = 4;
//        ProteinL = [|"name-1-4"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|-1.495405434; -0.1714587399; 0.01373039242; -0.4199918722; 0.5844246814; 1.488700972|];};
//    {ID = 5;
//        ProteinL = [|"name-2-5"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.170334643; 1.092401338; -0.004731922306; -0.9101257786; -0.08223106729; -1.265647213|];};
//    {ID = 6;
//        ProteinL = [|"name-2-6"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.752335164; 0.2229842076; 0.02571703711; -0.06511967118; -0.9478241777; -0.9880925599|];};
//    {ID = 7;
//        ProteinL = [|"name-2-7"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.241509293; 0.9578576745; 0.354962629; -0.8194239235; -0.5219383865; -1.212967286|];};
//    {ID = 8;
//        ProteinL = [|"name-2-8"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.582923146; 0.8160805522; -0.3156114193; -0.4361005577; -0.5221598957; -1.125131826|];};
//    {ID = 9;
//        ProteinL = [|"name-2-9"|];
//        OriginalBin = [|"root"|];
//        BinL = [|"root"|];
//        dataL =
//        [|1.25210773; 1.086679189; -0.2835872432; -1.133427748; 0.01944862816; -0.9412205564|];}|]

//let matrix1 = data1 |> distMatrixWeightedOf distanceMatrixWeighted None
//let matrix2 = data2 |> distMatrixWeightedOf distanceMatrixWeighted None
//let matrix3 = data3 |> distMatrixWeightedOf distanceMatrixWeighted None
//let matrix4 = data4 |> distMatrixWeightedOf distanceMatrixWeighted None

//drawKinetik data1 [|1.; 2.; 3.; 4.; 5.; 6.|] "1. dataset: low noise, 5-5" |> Chart.Show
//drawKinetik data2 [|1.; 2.; 3.; 4.; 5.; 6.|] "2. dataset: moderate noise, 5-5" |> Chart.Show
//drawKinetik data3 [|1.; 2.; 3.; 4.; 5.; 6.|] "3. dataset: moderate noise, 2-2-5" |> Chart.Show
//drawKinetik data4 [|1.; 2.; 3.; 4.; 5.; 6.|] "4. dataset: moderate noise, gain is better than clustering" |> Chart.Show

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//let ssn1 = applySSNcombi 10 data1
//let ssn2 = applySSNcombi 10 data2
//let ssn3 = applySSNcombi 10 data3
//let ssn4 = applySSNcombi 10 data4

////// direct kMean clustering with KKZ and all possible k [2 .. 10] and determine the best configuration by Gain 

//let kmeanGroupsKKZ (k: int) (children: Map<string,Types.Item []> ) =

//    let data = children |> Map.toArray |> Array.map (fun (s,ar) -> (s,ar |> Array.map (fun p -> p.dataL) |> MatrixTopLevelOperators.matrix))

//    let clusters = ML.Unsupervised.IterativeClustering.compute PreClusterFunctions.distMatrix (PreClusterFunctions.initCgroupsKKZ) PreClusterFunctions.updateCentroid data k

//    data
//    |> Array.map (fun list -> (clusters.Classifier list |> fst),list )
//    |> Array.groupBy (fst)
//    |> Array.map (fun (_,list) -> 
//        list 
//        |> Array.map (fun (_,(bin,_)) -> (bin, children |> Map.find bin))
//        )

//let gainGet matrix f clusters  = 
//    ClusterCheck.checkSgain f matrix clusters

//let kmeanCollectedFunction f data matrix = 
//    let children = 
//        data
//        |> Array.map (fun i -> (i.ProteinL.[0],[|i|]))
//        |> Map.ofArray
//    let clusterSet =
//        [|2 .. children.Count|]
//        |> Array.map (fun k -> 
//            let clusters =
//                kmeanGroupsKKZ k children
//            let gain = clusters |> Array.map (Array.map (fun (s,c) -> c) >> Array.concat) |> gainGet matrix f
//            (k, gain, clusters))
//    clusterSet
//    |> Array.maxBy (fun (k,g,c) -> g)

//let kmeanCollection f data matrix = 
//    let children = 
//        data
//        |> Array.map (fun i -> (i.ProteinL.[0],[|i|]))
//        |> Map.ofArray
//    let clusterSet =
//        [|2 .. children.Count|]
//        |> Array.map (fun k -> 
//            let clusters =
//                kmeanGroupsKKZ k children
//            let gain = clusters |> Array.map (Array.map (fun (s,c) -> c) >> Array.concat) |> gainGet matrix f
//            (k, gain, clusters))
//    clusterSet


/////

//let f = (SSN.getStepGainNodeSetnR 10)

//let clustering1 = kmeanCollectedFunction f data1 matrix1
//let clustering2 = kmeanCollectedFunction f data2 matrix2
//let clustering3 = kmeanCollectedFunction f data3 matrix3
//let clustering4 = kmeanCollectedFunction f data4 matrix4

//let clusteringSet1 = kmeanCollection f data1 matrix1
//let clusteringSet2 = kmeanCollection f data2 matrix2
//let clusteringSet3 = kmeanCollection f data3 matrix3
//let clusteringSet4 = kmeanCollection f data4 matrix4

//// no difference, easy
//ssn1.GroupGain                      //41.41360694
//clustering1 |> (fun (_,b,_) -> b)   //41.41360694

//// Different! 3 clusters against 6
//ssn2.GroupGain                      //25.33327168
//clustering2 |> (fun (_,b,_) -> b)   //24.09141619

//// no difference, easy
//ssn3.GroupGain                      //26.02960987
//clustering3 |> (fun (_,b,_) -> b)   //26.02960987

//// Different! 3 clusters against 2
//ssn4.GroupGain                      // 29.59737416
//clustering4 |> (fun (_,b,_) -> b)   // 28.68860158

//////////////////////////////////////////////////////////////////////// first check different gain formulas, with simplified IC

//////// IC = nCurrent / nRoot
//let f_ setNR dCurrSum dPredSum numberCurr numberRoot =
//    let nC = float numberCurr
//    let nR = float setNR
//    let deltaDist = dPredSum - dCurrSum
//    let deltaSpec = nC/nR //-((nC/nR)*log2(nC/nR)+((nR-nC)/nR)*log2((nR-nC)/nR))
//    if numberCurr=numberRoot then
//        0.
//    else
//        deltaDist*deltaSpec

//let ssn1_ = Functions.SSN.createTree (f_ 10) (None) Types.Mode.SSN_combi data1 
//let ssn2_ = Functions.SSN.createTree (f_ 10) (None) Types.Mode.SSN_combi data2
//let ssn3_ = Functions.SSN.createTree (f_ 10) (None) Types.Mode.SSN_combi data3
//let ssn4_ = Functions.SSN.createTree (f_ 10) (None) Types.Mode.SSN_combi data4


//let clustering1_ = kmeanCollectedFunction (f_ 10) data1 matrix1
//let clustering2_ = kmeanCollectedFunction (f_ 10) data2 matrix2
//let clustering3_ = kmeanCollectedFunction (f_ 10) data3 matrix3
//let clustering4_ = kmeanCollectedFunction (f_ 10) data4 matrix4


//// no difference
//ssn1_.GroupGain                      // 20.70680347
//clustering1_ |> (fun (_,b,_) -> b)   // 20.70680347

//// no difference
//ssn2_.GroupGain                      // 11.95730345
//clustering2_ |> (fun (_,b,_) -> b)   // 11.95730345

//// no difference, but both have only 2 clusters!            => bad IC choice
//ssn3_.GroupGain                      // 12.44000862
//clustering3_ |> (fun (_,b,_) -> b)   // 12.44000862

//// no difference
//ssn4_.GroupGain                      // 14.34430079
//clustering4_ |> (fun (_,b,_) -> b)   // 14.34430079

//////// IC = -(nC/nR)*log2(nC/nR)
//let f__ setNR dCurrSum dPredSum numberCurr numberRoot =
//    let nC = float numberCurr
//    let nR = float setNR
//    let deltaDist = dPredSum - dCurrSum
//    let deltaSpec = -(nC/nR)*log2(nC/nR)//+((nR-nC)/nR)*log2((nR-nC)/nR))
//    if numberCurr=numberRoot then
//        0.
//    else
//        deltaDist*deltaSpec

//let ssn1__ = Functions.SSN.createTree (f__ 10) (None) Types.Mode.SSN_combi data1 
//let ssn2__ = Functions.SSN.createTree (f__ 10) (None) Types.Mode.SSN_combi data2
//let ssn3__ = Functions.SSN.createTree (f__ 10) (None) Types.Mode.SSN_combi data3
//let ssn4__ = Functions.SSN.createTree (f__ 10) (None) Types.Mode.SSN_combi data4

//let clustering1__ = kmeanCollectedFunction (f__ 10) data1 matrix1
//let clustering2__ = kmeanCollectedFunction (f__ 10) data2 matrix2
//let clustering3__ = kmeanCollectedFunction (f__ 10) data3 matrix3
//let clustering4__ = kmeanCollectedFunction (f__ 10) data4 matrix4

//// ssn has 4 clusters!                                 => bad IC choice
//ssn1__.GroupGain                      // 21.24604748
//clustering1__ |> (fun (_,b,_) -> b)   // 20.70680347

//// no difference, but both have 6 clusters!            => bad IC choice
//ssn2__.GroupGain                      // 15.65092915
//clustering2__ |> (fun (_,b,_) -> b)   // 15.65092915

//// no difference, but both have 5 clusters!            => bad IC choice
//ssn3__.GroupGain                      // 15.70557879
//clustering3__ |> (fun (_,b,_) -> b)   // 15.70557879

//// no difference
//ssn4__.GroupGain                      // 17.1667973
//clustering4__ |> (fun (_,b,_) -> b)   // 17.1667973

//////// IC = nC/(nR/2. + nC) (michaelis-menten eq. with half root size as a Km)
//let f___ setNR dCurrSum dPredSum numberCurr numberRoot =
//    let nC = float numberCurr
//    let nR = float setNR
//    let deltaDist = dPredSum - dCurrSum
//    let deltaSpec = nC/(nR/2. + nC)//+((nR-nC)/nR)*log2((nR-nC)/nR))
//    if numberCurr=numberRoot then
//        0.
//    else
//        deltaDist*deltaSpec

//let ssn1___ = Functions.SSN.createTree (f___ 10) (None) Types.Mode.SSN_combi data1 
//let ssn2___ = Functions.SSN.createTree (f___ 10) (None) Types.Mode.SSN_combi data2
//let ssn3___ = Functions.SSN.createTree (f___ 10) (None) Types.Mode.SSN_combi data3
//let ssn4___ = Functions.SSN.createTree (f___ 10) (None) Types.Mode.SSN_combi data4

//let clustering1___ = kmeanCollectedFunction (f___ 10) data1 matrix1
//let clustering2___ = kmeanCollectedFunction (f___ 10) data2 matrix2
//let clustering3___ = kmeanCollectedFunction (f___ 10) data3 matrix3
//let clustering4___ = kmeanCollectedFunction (f___ 10) data4 matrix4

//// no difference,                              
//ssn1___.GroupGain                      // 20.70680347
//clustering1___ |> (fun (_,b,_) -> b)   // 20.70680347

//// no difference,         
//ssn2___.GroupGain                      // 11.95730345
//clustering2___ |> (fun (_,b,_) -> b)   // 11.95730345

//// no difference, but both have 4 clusters!            => bad IC choice
//ssn3___.GroupGain                      // 12.26356774
//clustering3___ |> (fun (_,b,_) -> b)   // 12.26356774

//// no difference
//ssn4___.GroupGain                      // 14.34430079
//clustering4___ |> (fun (_,b,_) -> b)   // 14.34430079

