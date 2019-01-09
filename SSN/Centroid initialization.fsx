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
//open TestData
//open GePhi
//open Types
//open Auxilliary
//open Plots
//open FunctionsExp

//#time

/////// Convert centroids into an initial scheme for K-Mean-Swap
////let centroidsToScheme (input: Item [] []) (centroid: matrix []) (scheme: int []) : ((string*(Item [])) []) =

////    let rec loop itemsRest groupID =
////        [|if groupID=scheme.Length then   
////            let binName = itemsRest |> Array.concat |> Array.map (fun p -> p.ID) |> Array.sort |> fun i -> String.Join("|",i)
////            yield (binName,itemsRest |> Array.concat)
////        else    
////            let (itemsCluster,itemsNew) = 
////                itemsRest 
////                |> Array.sortByDescending 
////                    (fun i -> PreClusterFunctions.pairwiseCorrAverage centroid.[groupID] (i |> Array.map (fun x -> x.dataL) |> MatrixTopLevelOperators.matrix))
////                |> Array.splitAt scheme.[groupID]
////            let binName = itemsCluster |> Array.concat |> Array.map (fun p -> p.ID) |> Array.sort |> fun i -> String.Join("|",i)
////            yield (binName,itemsCluster |> Array.concat)
////            yield! loop itemsNew (groupID+1)
////        |]
////    loop (input) 0

////let clustersToCentroidMatrix (clusters: Item [] [] []) : matrix [] =
////    clusters
////    |> Array.map (fun cluster -> 
////        cluster
////        |> Array.map (fun group ->
////            group |> Array.map (fun x -> x.dataL)
////            )
////        |> Array.concat
////        |> MatrixTopLevelOperators.matrix
////        )

//let centroidsToMatrix (centroids: Item [] list) : matrix [] =


//    centroids
//    |> List.toArray
//    |> Array.map (fun centroid ->
//        centroid |> Array.map (fun x -> x.dataL) |> MatrixTopLevelOperators.matrix
//        )


//let centroidsToClusterSigns (input: Item [] []) (centroids': Item [] list) (scheme: int []) : int [] =
//    let centroids = centroidsToMatrix centroids'
//    let rec loop itemsRest groupID =
//        [|if groupID=scheme.Length then   
//            yield (itemsRest)// |> Array.concat)
//        else    
//            let (itemsCluster,itemsNew) = 
//                itemsRest 
//                |> Array.sortByDescending 
//                    (fun i -> PreClusterFunctions.pairwiseCorrAverage centroids.[groupID] (i |> Array.map (fun x -> x.dataL) |> MatrixTopLevelOperators.matrix))
//                |> Array.splitAt scheme.[groupID]
            
//            yield (itemsCluster)// |> Array.concat)
//            yield! loop itemsNew (groupID+1)
//        |]
//    loop (input) 0
//    |> Array.concat
//    |> Array.map (fun x -> Array.FindIndex(input, fun i -> i=x))


////let kkzSingle (data: Item []) k =
////    let centroid1 =
////        data |> Array.maxBy (fun x -> sqrt ( x.dataL |> Array.sumBy (fun i -> i*i)))
////    let LeaveData d c =
////        d |> Array.removeIndex (Array.FindIndex<Item>(d, fun x -> x=c))
////    let rec loop dataRest kRest centroids =
////        if kRest=1 then   
////            centroids
////        else    
////            let newC = 
////                dataRest 
////                |> Array.map  (fun p -> p, centroids |> List.map (fun c -> SSN.weightedEuclidean None p.dataL c.dataL) |> List.min )
////                |> Array.maxBy snd 
////                |> fst
////            loop (LeaveData dataRest newC) (kRest-1) (newC::centroids)
////    loop (LeaveData data centroid1) k [centroid1]

///// give a list of centroids as k the most distant elements of the dataset     
//let kkz (data: Item [] []) k =
//    let centroid1 =
//        data 
//        |> Array.maxBy 
//            (fun x -> 
//                x 
//                |> Array.map (fun x -> x.dataL) 
//                |> JaggedArray.transpose 
//                |> Array.map (Array.average) 
//                |> fun xx ->  sqrt ( xx |> Array.sumBy (fun i -> i*i))
//            )
//    let LeaveData d c =
//        d |> Array.removeIndex (Array.FindIndex<Item []>(d, fun x -> x=c))
//    let toMatrix =
//        Array.map (fun i -> i.dataL)
//        >> MatrixTopLevelOperators.matrix

//    let rec loop dataRest kRest centroids =
//        if kRest=1 then   
//            centroids
//        else    
//            let newC = 
//                dataRest 
//                |> Array.map (fun p -> 
//                    p, centroids |> List.map (fun c -> PreClusterFunctions.pairwiseCorrAverage (toMatrix p)  (toMatrix c)) |> List.min )
//                |> Array.maxBy snd 
//                |> fst
//            loop (LeaveData dataRest newC) (kRest-1) (newC::centroids)
//    loop (LeaveData data centroid1) k [centroid1]

////let varPart_Single (data: Item []) k =
////    let sse (cluster: Item []) =
////        let centroid = cluster |> Array.map (fun x -> x.dataL) |> MatrixTopLevelOperators.matrix
////        let dist (a: Item) (b: matrix) =
////            PreClusterFunctions.pairwiseCorrAverage ([a.dataL] |> MatrixTopLevelOperators.matrix) b
////        cluster
////        |> Array.map 
////            (fun i -> (dist i centroid)*(dist i centroid))
////        |> Array.average
////    let split (cluster: Item []) =
////        let featureN = 
////            cluster 
////            |> Array.map (fun x -> x.dataL) 
////            |> MatrixTopLevelOperators.matrix 
////            |> Matrix.Generic.enumerateColumnWise Seq.var
////            |> Seq.mapi (fun id x -> (id,x))
////            |> Seq.maxBy snd
////            |> fst
////        let featureMean =
////            cluster
////            |> Array.map (fun x -> x.dataL.[featureN])
////            |> Array.average
////        cluster
////        |> Array.partition (fun x -> x.dataL.[featureN]>featureMean)
////    let pq = IndexPriorityQueue<float>(k*2-1) // put a cluster there or SSE of a cluster or SSE*cluster????
////    let clusters: Item [] [] = Array.create (k*2-1) [||]
////    pq.Insert 0 (sse data)
////    clusters.[0] <- data
////    [|1 .. (k-1)|] 
////    |> Array.iter (fun ik -> 
////        let loosest = clusters.[pq.HeapItemIndex 1]
////        let newCl = split loosest
////        pq.Pop() |> ignore
////        pq.Insert (2*ik) (sse (fst newCl))
////        pq.Insert (2*ik-1) (sse (snd newCl))
////        clusters.[2*ik] <- (fst newCl)
////        clusters.[2*ik-1] <- (snd newCl)
////        )
////    [|1 .. k|]
////    |> Array.map (fun x -> clusters.[pq.HeapItemIndex x])


////let varPart (data: Item [] []) k =
////    let sse (cluster: Item [] []) =
////        let centroid = cluster |> Array.map (fun x -> x |> Array.map (fun i -> i.dataL)) |> Array.concat |> MatrixTopLevelOperators.matrix
////        let dist (a: Item []) (b: matrix) =
////            PreClusterFunctions.pairwiseCorrAverage (a |> Array.map (fun i -> i.dataL) |> MatrixTopLevelOperators.matrix) b
////        cluster
////        |> Array.map 
////            (fun i -> (dist i centroid)*(dist i centroid))
////        |> Array.average
////    let split (cluster: Item [] []) =
////        let featureN = 
////            cluster 
////            |> Array.map (fun x -> x |> Array.map (fun i -> i.dataL) |> JaggedArray.transpose |> Array.map Array.average )
////            |> MatrixTopLevelOperators.matrix 
////            |> Matrix.Generic.enumerateColumnWise Seq.var
////            |> Seq.mapi (fun id x -> (id,x))
////            |> Seq.maxBy snd
////            |> fst
////        let featureMean =
////            cluster
////            |> Array.map (fun x -> x |> Array.map (fun i -> i.dataL.[featureN]) |> Array.average)
////            |> Array.average
////        cluster
////        |> Array.partition (fun x -> (x |> Array.map (fun i -> i.dataL.[featureN]) |> Array.average)>featureMean)
////    let pq = IndexPriorityQueue<float>(k*2-1) // put a cluster there or SSE of a cluster or SSE*cluster????
////    let clusters: Item [] [] [] = Array.create (k*2-1) [||]
////    pq.Insert 0 (sse data)
////    clusters.[0] <- data
////    [|1 .. (k-1)|] 
////    |> Array.iter (fun ik -> 
////        let loosest = clusters.[pq.HeapItemIndex 1]
////        let newCl = split loosest
////        pq.Pop() |> ignore
////        pq.Insert (2*ik) (sse (fst newCl))
////        pq.Insert (2*ik-1) (sse (snd newCl))
////        clusters.[2*ik] <- (fst newCl)
////        clusters.[2*ik-1] <- (snd newCl)
////        )
////    [|1 .. k|]
////    |> Array.map (fun x -> clusters.[pq.HeapItemIndex x])


//// testing

////let data6t_35 = 
////    ChlamyProteome.dataAll
////    |> Array.filter (fun x -> x.BinL.[0]="35")
////    |> Array.mapi (fun id x -> {x with ID=id})
////    |> Array.splitAt 10
////    |> fst

//////varPart_Single data6t_35 3
////kkz (data6t_35 |> Array.map (fun i -> [|i|])) 3
//////varPart (data6t_35 |> Array.map (fun i -> [|i|])) 3

//// how we have to incorporate a centroid initialization in K-Mean and centroid initialization converter in K-Mean-Swap 
//// maybe omit k-mean at all in favour of Hierarchical clustering

//let intiCgroupsKKZ (data: (string*matrix) array) k =
//    let centroid1 =
//        data 
//        |> Array.maxBy 
//            (fun x -> 
//                x 
//                |> snd
//                |> Matrix.enumerateColumnWise (Seq.mean) 
//                |> fun xx ->  sqrt ( xx |> Seq.sumBy (fun i -> i*i))
//            )
//    let LeaveData d c =
//        d |> Array.removeIndex (Array.FindIndex<string*matrix>(d, fun x -> x=c))

//    let rec loop dataRest kRest centroids =
//        if kRest=1 then   
//            centroids
//        else    
//            let newC = 
//                dataRest 
//                |> Array.map (fun (s,p) -> 
//                    (s,p), centroids |> List.map (fun (sc,c) -> PreClusterFunctions.pairwiseCorrAverage (p)  (c)) |> List.min )
//                |> Array.maxBy snd 
//                |> fst
//            loop (LeaveData dataRest newC) (kRest-1) (newC::centroids)

//    loop (LeaveData data centroid1) k [centroid1]
//    |> List.toArray
    
//let kmeanGroupsKKZ (k: int) weight (children: Map<string,Types.Item []> ) : Map<string,Types.Item []> =

//        let data = children |> Map.toArray |> Array.map (fun (s,ar) -> (s,ar |> Array.map (fun p -> p.dataL) |> MatrixTopLevelOperators.matrix))

//        let clusters = 
//            let c1 = ML.Unsupervised.IterativeClustering.compute FunctionsExp.PreClusterFunctions.distMatrix (intiCgroupsKKZ) FunctionsExp.PreClusterFunctions.updateCentroid data k
//            let x1 = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c1
//            [|1 .. 20|]
//            |> Array.fold (fun (disp,best) x -> 
//                let c = ML.Unsupervised.IterativeClustering.compute FunctionsExp.PreClusterFunctions.distMatrix (intiCgroupsKKZ) FunctionsExp.PreClusterFunctions.updateCentroid data k
//                let x = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c
//                if x<disp then
//                    (x,c)
//                else
//                    (disp,best) ) (x1,c1)
//            |> snd

//        data
//        |> Array.map (fun list -> (clusters.Classifier list |> fst),list )
//        |> Array.groupBy (fst)
//        |> Array.map (fun (cID,list) ->
//            let binName = list |> Array.map (fun (cID,(bin,p)) -> bin) |> Array.sort |> fun i -> String.Join("|",i) 
//            let items = list |> Array.map (fun (cID,(bin,p)) -> children |> Map.find bin) |> Array.concat
//            (binName,items))
//        |> Map.ofArray


////let intiCgroups_VarPart (data: (string*matrix) array) k =
////    let sse (cluster: (string*matrix) []) =
////        let centroid = cluster |> Array.map (snd >> Matrix.enumerateRowWise (fun x -> x |> Seq.toArray) >> Seq.toArray) |> Array.concat |> MatrixTopLevelOperators.matrix
////        let dist (a: matrix) (b: matrix) =
////            PreClusterFunctions.pairwiseCorrAverage (a) b
////        cluster
////        |> Array.map 
////            (fun i -> (dist (snd i) centroid)*(dist (snd i) centroid))
////        |> Array.average
////    let split (cluster: (string*matrix) []) =
////        let featureN = 
////            cluster 
////            |> Array.map (snd >> Matrix.enumerateColumnWise (fun x -> x |> Seq.mean) )
////            |> MatrixTopLevelOperators.matrix 
////            |> Matrix.Generic.enumerateColumnWise Seq.var
////            |> Seq.mapi (fun id x -> (id,x))
////            |> Seq.maxBy snd
////            |> fst
////        let featureMean =
////            cluster
////            |> Array.map (fun (s,x) -> Matrix.getCol x featureN |> Vector.mean)
////            |> Array.average
////        cluster
////        |> Array.partition (fun (s,x) -> ( Matrix.getCol x featureN |> Vector.mean ) > featureMean)
////    let pq = IndexPriorityQueue<float>(k*2-1) // put a cluster there or SSE of a cluster or SSE*cluster????
////    let clusters: (string*matrix) [] [] = Array.create (k*2-1) [||]
////    pq.Insert 0 (sse data)
////    clusters.[0] <- data
////    [|1 .. (k-1)|] 
////    |> Array.iter (fun ik -> 
////        let loosest = clusters.[pq.HeapItemIndex 1]
////        let newCl = split loosest
////        pq.Pop() |> ignore
////        pq.Insert (2*ik) (sse (fst newCl))
////        pq.Insert (2*ik-1) (sse (snd newCl))
////        clusters.[2*ik] <- (fst newCl)
////        clusters.[2*ik-1] <- (snd newCl)
////        )
////    [|1 .. k|]
////    |> Array.map (fun x -> clusters.[pq.HeapItemIndex x])
////    |> Array.map (fun cluster -> 
////        let bin = String.Join("|",(cluster |> Array.map fst))
////        let clusterMatrix = cluster |> Array.map (snd >> Matrix.enumerateRowWise (fun x -> x |> Seq.toArray) >> Seq.toArray) |> Array.concat |> MatrixTopLevelOperators.matrix
////        (bin,clusterMatrix))

////let kmeanGroupsVP (k: int) weight (children: Map<string,Types.Item []> ) : Map<string,Types.Item []> =

////        let data = children |> Map.toArray |> Array.map (fun (s,ar) -> (s,ar |> Array.map (fun p -> p.dataL) |> MatrixTopLevelOperators.matrix))

////        let clusters = 
////            let c1 = ML.Unsupervised.IterativeClustering.compute FunctionsExp.PreClusterFunctions.distMatrix (intiCgroups_VarPart) FunctionsExp.PreClusterFunctions.updateCentroid data k
////            let x1 = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c1
////            [|1 .. 20|]
////            |> Array.fold (fun (disp,best) x -> 
////                let c = ML.Unsupervised.IterativeClustering.compute FunctionsExp.PreClusterFunctions.distMatrix (intiCgroups_VarPart) FunctionsExp.PreClusterFunctions.updateCentroid data k
////                let x = ML.Unsupervised.IterativeClustering.DispersionOfClusterResult c
////                if x<disp then
////                    (x,c)
////                else
////                    (disp,best) ) (x1,c1)
////            |> snd

////        data
////        |> Array.map (fun list -> (clusters.Classifier list |> fst),list )
////        |> Array.groupBy (fst)
////        |> Array.map (fun (cID,list) ->
////            let binName = list |> Array.map (fun (cID,(bin,p)) -> bin) |> Array.sort |> fun i -> String.Join("|",i) 
////            let items = list |> Array.map (fun (cID,(bin,p)) -> children |> Map.find bin) |> Array.concat
////            (binName,items))
////        |> Map.ofArray


///// testing


////let data35 = 
////    ChlamyProteome.dataAll
////    |> Array.filter (fun x -> x.BinL.[0]="35")
////    |> Array.mapi (fun id x -> {x with ID=id})

////let data35toPre = 
////    ChlamyProteome.dataAll
////    |> Array.filter (fun x -> x.BinL.[0]="35" && x.BinL.[1]="2")
////    |> Array.mapi (fun id x -> {x with ID=id})

////let data35toPreMap = 
////    data35toPre
////    |> Array.map (fun x -> (string x.ID, [|x|]))
////    |> Map.ofArray

////let matrix = 
////    data35toPre
////    |> SSN.distMatrixWeightedOf SSN.distanceMatrixWeighted None

////let gainFn = SSN.getStepGainNodeSetnR 300


////let newClustersKKZ = kmeanGroupsKKZ 19 None data35toPreMap

////let preG_KKZ = 
////    newClustersKKZ 
////    |> Map.toArray 
////    |> Array.sumBy (fun (s,items) -> SSN.getStepGainFn gainFn (items |> SSN.groupIDFn) (data35toPre |> SSN.groupIDFn) data35toPre.Length matrix)
////printfn "pre-Cluster G = %f" preG_KKZ

////let applySSN_Reduction_KKZ_random data setN = SSN.createTreeReducing 19 (SSN.getStepGainNodeSetnR setN) (KMeanSwapFunctions.kmeanSwapShuffle setN 1) kmeanGroupsKKZ None SSN.SSN data

////applySSN_Reduction_KKZ_random data35 300

////let f_KKZ =
////    [for a in [1 .. 50] ->    
////        let x = applySSN_Reduction_KKZ_random data35 300
////        printfn "G=%f" x.GroupGain
////        x
////        ]

////let best_KKZ = f_KKZ |> List.maxBy (fun i -> i.GroupGain)
////best_KKZ.GroupGain


////let newClustersVP = kmeanGroupsVP 30 None data35toPreMap

////let preG_VP = 
////    newClustersVP
////    |> Map.toArray 
////    |> Array.sumBy (fun (s,items) -> SSN.getStepGainFn gainFn (items |> SSN.groupIDFn) (data35toPre |> SSN.groupIDFn) data35toPre.Length matrix)
////printfn "pre-Cluster G = %f" preG_VP

////let applySSN_Reduction_VP_random data setN = SSN.createTreeReducing 19 (SSN.getStepGainNodeSetnR setN) (KMeanSwapFunctions.kmeanSwapShuffle setN 1) kmeanGroupsVP None SSN.SSN data

////applySSN_Reduction_VP_random data35 300

////let f_VP =
////    [for a in [1 .. 50] ->    
////        let x = applySSN_Reduction_VP_random data35 300
////        printfn "G=%f" x.GroupGain
////        x
////        ]

////let best_VP = f_VP |> List.maxBy (fun i -> i.GroupGain)
////best_VP.GroupGain

//let f setNR dCurrSum dPredSum numberCurr numberRoot =
//        let nC = float numberCurr
//        let nR = float numberRoot
//        let deltaDist = dPredSum - dCurrSum
//        let deltaSpec = nC/nR // nC/(nR/2.+nC)//-((nC/nR)*log2(nC/nR)+((nR-nC)/nR)*log2((nR-nC)/nR))
//        if numberCurr=numberRoot then
//            0.
//        else
//            deltaDist*deltaSpec

///// KMean Swap with deterministic seeding, with criteria based on Gain formula
//let kMeanSwapDet nSet f matrixSingletons depth (initialGrouping': Map<string,(Types.Item array)>) : (Map<string,Types.Item array>) [] =
//    let initialGrouping = initialGrouping' |> Map.toArray
//    let l = initialGrouping.Length

//    let parentGroup = initialGrouping |> Array.map (fun (k,p) -> p) |> Array.concat

//    let mutable iterN : int = 0

//    let schemes = seq [2..l] |> Seq.map (SSN.schemeGenerator l) |> Seq.concat
//    let centroidsSet = [2..l] |> List.map (kkz (initialGrouping |> Array.map snd))

//    let gainDiff listA listB : float =
//        let gFn (current: Types.Item array) = 
//            SSN.getStepGainFn f (SSN.groupIDFn current) (SSN.groupIDFn parentGroup) parentGroup.Length matrixSingletons
//        (gFn (Array.append listA listB)) //- (gFn listA) - (gFn listB)

//    let matrix =     
//        let data =
//            initialGrouping
//            |> Array.map snd
//        let m = Array2D.zeroCreate (data.Length) (data.Length)
//        for rowI in 0..data.Length-1 do
//            for colI in 0..rowI do
//                let tmp = if rowI=colI then nan else gainDiff data.[rowI] data.[colI]
//                m.[colI,rowI] <- tmp
//                m.[rowI,colI] <- tmp
//        m

//    let processInt (intListList: int [] []) =
//        intListList
//        |> Array.map (fun i ->
//                                if i.Length=1 then
//                                    initialGrouping.[i.[0]]
//                                else
//                                    Array.fold (fun (key,gr) ii ->
//                                    let (binKey, value) = initialGrouping.[ii]
//                                    ((sprintf "%s|%s" key binKey), Array.append value gr )
//                                    ) ("",[||]) i)
//        |> Array.map (fun (newBin,protA) ->
//            (newBin, protA |> Array.map (fun prot ->
//                    {prot with BinL = Array.append prot.BinL.[0..(depth)] [|newBin|]}))) // updating the bin labels            
//        |> Map.ofArray

//    let rec swap (matrix: float [,]) (clusterN: int) (currentConf: int [] []) n downIterLast bestSoFar (lastGain: KMeanSwapFunctions.GainComponents [] []) (swappedGroups: int []) (swappedElements: int []) = /// K-Mean-Swap for matrices
            
//        let newGain =
//            if swappedGroups=[||] then
//                lastGain
//            else
//                lastGain
//                |> Array.mapi (fun idGroup group ->
//                    if (idGroup = swappedGroups.[0]) then
//                        group
//                        |> Array.mapi (fun idElement element ->

//                            let currentGroupIDs = currentConf.[idGroup] |> Array.map (fun i -> initialGrouping.[i] |> snd) |> Array.concat  |> SSN.groupIDFn
//                            let currentElement = initialGrouping.[currentConf.[idGroup].[idElement]] |> snd |> SSN.groupIDFn

//                            if idElement = swappedElements.[0] then
//                                {ParentDiss = SSN.dSumFn currentElement (parentGroup |> SSN.groupIDFn) matrixSingletons;
//                                CurrentDiss = SSN.dSumFn currentElement currentGroupIDs matrixSingletons}
//                            else
//                                {element with CurrentDiss = SSN.dSumFn currentElement currentGroupIDs matrixSingletons})
//                    elif (idGroup = swappedGroups.[1]) then
//                        group
//                        |> Array.mapi (fun idElement element ->
                        
//                            let currentGroupIDs = currentConf.[idGroup] |> Array.map (fun i -> initialGrouping.[i] |> snd) |> Array.concat  |> SSN.groupIDFn
//                            let currentElement = initialGrouping.[currentConf.[idGroup].[idElement]] |> snd |> SSN.groupIDFn

//                            if idElement = swappedElements.[1] then
//                                {ParentDiss = SSN.dSumFn currentElement (parentGroup |> SSN.groupIDFn) matrixSingletons;
//                                    CurrentDiss = SSN.dSumFn currentElement currentGroupIDs matrixSingletons}
//                            else
//                                {element with CurrentDiss = SSN.dSumFn currentElement currentGroupIDs matrixSingletons})
//                    else
//                        group
//                    )

//        let newGainValue =
//            KMeanSwapFunctions.gainValueFn nSet newGain

//        let (downIter, best) =
//            if (newGainValue < (fst bestSoFar)) then
//                (downIterLast+1, bestSoFar)
//            else
//                (downIterLast, (newGainValue,currentConf))
        
//        let cluster = currentConf.[clusterN]

//        let outside =
//            currentConf
//            |> Array.concat
//            |> Array.except cluster
//        let (idOut, idIn, minDist) =
//            outside
//            |> Array.map (fun iOut -> Array.map (fun iIn -> (iOut, iIn, matrix.[iOut,iIn])) cluster)
//            |> Array.concat
//            |> Array.maxBy (fun (_, _, distance) -> distance) // find the closest item from outside,
//        let (farItem, _, maxDist) =
//            cluster
//            |> Array.filter (fun iIn -> iIn<>idIn)
//            |> Array.map (fun iIn -> (iIn, idIn, matrix.[iIn,idIn]))
//            |> Array.minBy (fun (_, _, distance) -> distance) // find the farthest to the inside item within a cluster
            
//        if (n >= 0) // size of the cluster as a limit for amount of iterations
//            && (downIter<2) // no more than 2 down-iterations pro cluster
//            && (maxDist) < (minDist) // there is a gain in swapping !!! try to figure out another criteria for switching
//            then
//                let idInsideItem = cluster |> Array.findIndex (fun id -> id = farItem)
//                let idOutsideCluster = currentConf |>  Array.findIndex (fun i -> i |> Array.exists (fun (ii) -> ii = idOut))
//                let idOutsideItem = Array.findIndex (fun id -> id = idOut) currentConf.[idOutsideCluster]
//                let temp = cluster.[idInsideItem]
//                currentConf.[clusterN].[idInsideItem] <- currentConf.[idOutsideCluster].[idOutsideItem]
//                currentConf.[idOutsideCluster].[idOutsideItem] <- temp
//                iterN <- iterN + 1
//                swap matrix clusterN currentConf (n-1) downIter best newGain [|clusterN;idOutsideCluster|] [|idInsideItem;idOutsideItem|]
//            else best

//    let fillSchemes (scheme: int []) : Map<string, Types.Item array> =
        
//        let initialSchemeFilling =
//            let centroids = centroidsSet.[scheme.Length-2]
//            centroidsToClusterSigns (initialGrouping |> Array.map snd) centroids scheme

//        let configurationIDs = // make a change here, to fill the scheme in order of the clustering result
//            scheme
//            |> Array.mapFold (fun before n -> ([|before .. (before+n-1)|],(before+n))) 0
//            |> fst
//            |> Array.map (fun i -> i |> Array.map (fun ii -> initialSchemeFilling.[ii]))
        
//        let result =
//            [|0 .. (scheme.Length-1)|]
//            |> Array.fold (fun (gainBest, bestSoFar) clusterN ->
//                if (scheme.[clusterN]>1) then
                    
//                    let gainComp =
//                        configurationIDs
//                        |> Array.map (fun idGroup ->
//                            Array.map (fun idElement ->
//                                let currentGroupIDs = idGroup |> Array.map (fun i -> initialGrouping.[i] |> snd) |> Array.concat  |> SSN.groupIDFn
//                                let currentElement = initialGrouping.[idElement] |> snd |> SSN.groupIDFn
//                                {KMeanSwapFunctions.ParentDiss = SSN.dSumFn currentElement (parentGroup |> SSN.groupIDFn) matrixSingletons;
//                                KMeanSwapFunctions.CurrentDiss = SSN.dSumFn currentElement currentGroupIDs matrixSingletons}) idGroup )

//                    swap matrix clusterN bestSoFar scheme.[clusterN] 0 (0.,[||]) gainComp [||] [||]
//                else
//                    (gainBest,bestSoFar)) (0., configurationIDs)
//            |> snd
            
//        result |> processInt
    
//    let r =
//        schemes
//        |> Seq.rev
//        |> List.ofSeq
//        |> List.map fillSchemes
//        |> Array.ofList

//    printfn "iterations done: %i" iterN

//    r


///// KMean Swap with deterministic seeding, with criteria based only on Distance to centroid
//let kMeanSwapDetCentr nSet f matrixSingletons depth (initialGrouping': Map<string,(Types.Item array)>) : (Map<string,Types.Item array>) [] =
//    let initialGrouping = initialGrouping' |> Map.toArray
//    let l = initialGrouping.Length

//    let parentGroup = initialGrouping |> Array.map (fun (k,p) -> p) |> Array.concat

//    let mutable iterN : int = 0

//    let schemes = seq [2..l] |> Seq.map (SSN.schemeGenerator l) |> Seq.concat
//    let centroidsSet = [2..l] |> List.map (kkz (initialGrouping |> Array.map snd))

//    let gainDiff listA listB : float =
//        //let gFn (current: Types.Item array) = 
//        //    SSN.getStepGainFn f (SSN.groupIDFn current) (SSN.groupIDFn parentGroup) parentGroup.Length matrixSingletons
//        //(gFn (Array.append listA listB)) //- (gFn listA) - (gFn listB)
//        let matr list =
//            list |> Array.map (fun i -> i.dataL) |> MatrixTopLevelOperators.matrix
//        - (PreClusterFunctions.pairwiseCorrMax (matr listA) (matr listB))

//    let matrix =     
//        let data =
//            initialGrouping
//            |> Array.map snd
//        let m = Array2D.zeroCreate (data.Length) (data.Length)
//        for rowI in 0..data.Length-1 do
//            for colI in 0..rowI do
//                let tmp = if rowI=colI then nan else gainDiff data.[rowI] data.[colI]
//                m.[colI,rowI] <- tmp
//                m.[rowI,colI] <- tmp
//        m

//    let processInt (intListList: int [] []) =
//        intListList
//        |> Array.map (fun i ->
//                                if i.Length=1 then
//                                    initialGrouping.[i.[0]]
//                                else
//                                    Array.fold (fun (key,gr) ii ->
//                                    let (binKey, value) = initialGrouping.[ii]
//                                    ((sprintf "%s|%s" key binKey), Array.append value gr )
//                                    ) ("",[||]) i)
//        |> Array.map (fun (newBin,protA) ->
//            (newBin, protA |> Array.map (fun prot ->
//                    {prot with BinL = Array.append prot.BinL.[0..(depth)] [|newBin|]}))) // updating the bin labels            
//        |> Map.ofArray

//    let rec swap (matrix: float [,]) (clusterN: int) (currentConf: int [] []) n downIterLast bestSoFar (lastGain: KMeanSwapFunctions.GainComponents [] []) (swappedGroups: int []) (swappedElements: int []) = /// K-Mean-Swap for matrices
            
//        let newGain =
//            if swappedGroups=[||] then
//                lastGain
//            else
//                lastGain
//                |> Array.mapi (fun idGroup group ->
//                    if (idGroup = swappedGroups.[0]) then
//                        group
//                        |> Array.mapi (fun idElement element ->

//                            let currentGroupIDs = currentConf.[idGroup] |> Array.map (fun i -> initialGrouping.[i] |> snd) |> Array.concat  |> SSN.groupIDFn
//                            let currentElement = initialGrouping.[currentConf.[idGroup].[idElement]] |> snd |> SSN.groupIDFn

//                            if idElement = swappedElements.[0] then
//                                {ParentDiss = SSN.dSumFn currentElement (parentGroup |> SSN.groupIDFn) matrixSingletons;
//                                CurrentDiss = SSN.dSumFn currentElement currentGroupIDs matrixSingletons}
//                            else
//                                {element with CurrentDiss = SSN.dSumFn currentElement currentGroupIDs matrixSingletons})
//                    elif (idGroup = swappedGroups.[1]) then
//                        group
//                        |> Array.mapi (fun idElement element ->
                        
//                            let currentGroupIDs = currentConf.[idGroup] |> Array.map (fun i -> initialGrouping.[i] |> snd) |> Array.concat  |> SSN.groupIDFn
//                            let currentElement = initialGrouping.[currentConf.[idGroup].[idElement]] |> snd |> SSN.groupIDFn

//                            if idElement = swappedElements.[1] then
//                                {ParentDiss = SSN.dSumFn currentElement (parentGroup |> SSN.groupIDFn) matrixSingletons;
//                                    CurrentDiss = SSN.dSumFn currentElement currentGroupIDs matrixSingletons}
//                            else
//                                {element with CurrentDiss = SSN.dSumFn currentElement currentGroupIDs matrixSingletons})
//                    else
//                        group
//                    )

//        let newGainValue =
//            KMeanSwapFunctions.gainValueFn nSet newGain

//        let (downIter, best) =
//            if (newGainValue < (fst bestSoFar)) then
//                (downIterLast+1, bestSoFar)
//            else
//                (downIterLast, (newGainValue,currentConf))
        
//        let cluster = currentConf.[clusterN]

//        let outside =
//            currentConf
//            |> Array.concat
//            |> Array.except cluster
//        let (idOut, idIn, minDist) =
//            outside
//            |> Array.map (fun iOut -> Array.map (fun iIn -> (iOut, iIn, matrix.[iOut,iIn])) cluster)
//            |> Array.concat
//            |> Array.maxBy (fun (_, _, distance) -> distance) // find the closest item from outside,
//        let (farItem, _, maxDist) =
//            cluster
//            |> Array.filter (fun iIn -> iIn<>idIn)
//            |> Array.map (fun iIn -> (iIn, idIn, matrix.[iIn,idIn]))
//            |> Array.minBy (fun (_, _, distance) -> distance) // find the farthest to the inside item within a cluster
            
//        if (n >= 0) // size of the cluster as a limit for amount of iterations
//            && (downIter<2) // no more than 2 down-iterations pro cluster
//            && (maxDist) < (minDist) // there is a gain in swapping !!! try to figure out another criteria for switching
//            then
//                let idInsideItem = cluster |> Array.findIndex (fun id -> id = farItem)
//                let idOutsideCluster = currentConf |>  Array.findIndex (fun i -> i |> Array.exists (fun (ii) -> ii = idOut))
//                let idOutsideItem = Array.findIndex (fun id -> id = idOut) currentConf.[idOutsideCluster]
//                let temp = cluster.[idInsideItem]
//                currentConf.[clusterN].[idInsideItem] <- currentConf.[idOutsideCluster].[idOutsideItem]
//                currentConf.[idOutsideCluster].[idOutsideItem] <- temp
//                iterN <- iterN + 1
//                swap matrix clusterN currentConf (n-1) downIter best newGain [|clusterN;idOutsideCluster|] [|idInsideItem;idOutsideItem|]
//            else best

//    let fillSchemes (scheme: int []) : Map<string, Types.Item array> =
        
//        let initialSchemeFilling =
//            let centroids = centroidsSet.[scheme.Length-2]
//            centroidsToClusterSigns (initialGrouping |> Array.map snd) centroids scheme

//        let configurationIDs = // make a change here, to fill the scheme in order of the clustering result
//            scheme
//            |> Array.mapFold (fun before n -> ([|before .. (before+n-1)|],(before+n))) 0
//            |> fst
//            |> Array.map (fun i -> i |> Array.map (fun ii -> initialSchemeFilling.[ii]))
        
//        let result =
//            [|0 .. (scheme.Length-1)|]
//            |> Array.fold (fun (gainBest, bestSoFar) clusterN ->
//                if (scheme.[clusterN]>1) then
                    
//                    let gainComp =
//                        configurationIDs
//                        |> Array.map (fun idGroup ->
//                            Array.map (fun idElement ->
//                                let currentGroupIDs = idGroup |> Array.map (fun i -> initialGrouping.[i] |> snd) |> Array.concat  |> SSN.groupIDFn
//                                let currentElement = initialGrouping.[idElement] |> snd |> SSN.groupIDFn
//                                {KMeanSwapFunctions.ParentDiss = SSN.dSumFn currentElement (parentGroup |> SSN.groupIDFn) matrixSingletons;
//                                KMeanSwapFunctions.CurrentDiss = SSN.dSumFn currentElement currentGroupIDs matrixSingletons}) idGroup )

//                    swap matrix clusterN bestSoFar scheme.[clusterN] 0 (0.,[||]) gainComp [||] [||]
//                else
//                    (gainBest,bestSoFar)) (0., configurationIDs)
//            |> snd
            
//        result |> processInt
    
//    let r =
//        schemes
//        |> Seq.rev
//        |> List.ofSeq
//        |> List.map fillSchemes
//        |> Array.ofList

//    printfn "iterations done: %i" iterN

//    r

////let kmeanSwap_Det setN power f matrixItems depth (map : Map<string,Types.Item array>) =
////    let items = map |> Map.toArray
////    [|for a in [1 .. power] do yield
////                                    items
////                                    |> Array.sortBy (fun (s,_) -> s)
////                                    |> (fun x ->
////                                            printfn "randomWalk %i" a
////                                            kMeanSwapDet setN f matrixItems depth x )
////                                    |> Seq.toArray|]
////    |> Array.stackVertical
////    |> Array2D.toJaggedArray
////    |> Array.map (fun a -> a |> Array.maxBy (fun m -> m |> Map.toArray |> Array.map (snd >> SSN.groupIDFn) |> KMeanSwapFunctions.gainstepcheck f matrixItems))

    
////let applySSN_Reduction_KKZ data setN pre_k = SSN.createTreeReducing pre_k (SSN.getStepGainNodeSetnR setN) (kmeanSwap_Det setN 1) kmeanGroupsKKZ None SSN.SSN data

////let best_KKZ_det = applySSN_Reduction_KKZ data35 300 26

////best_KKZ_det |> Tree.filterLeaves |> Array.map (fun i -> i.Length)

////let gvsprek = [2 .. 40] |> List.map (fun prek -> 
////    let tree = applySSN_Reduction_KKZ data35 300 prek
////    (prek, tree))

////GePhi.sendToGephiFromTreeParam (gvsprek.[27-2] |> snd)
 
////(gvsprek.[27-2] |> snd) |> Tree.filterLeaves |> Array.map (fun i -> i.Length)

////Chart.Line (gvsprek |> List.map (fun (x,y) -> x,y.GroupGain)) |> Chart.withX_AxisStyle "pre_k" |> Chart.withY_AxisStyle "root GroupGain" |> Chart.Show
////Chart.Line (gvsprek |> List.map (fun (x,z) -> x,(z.Children.Item "2").Children.Count)) |> Chart.withX_AxisStyle "pre_k" |> Chart.withY_AxisStyle "children number" |> Chart.Show

////[
////    Chart.Line((gvsprek |> List.map (fun (x,y) -> x,y.GroupGain)), ShowMarkers=true ,Name="root G_gain")
////    |> Chart.withAxisAnchor(Y=1);
////    Chart.Line((gvsprek |> List.map (fun (x,z) -> x,(z.Children.Item "2").Children.Count)), ShowMarkers=true , Name="children N")
////    |> Chart.withAxisAnchor(Y=2);
////]
////|> Chart.Combine
////|> Chart.withY_AxisStyle("root G_gain", Side=StyleParam.Side.Left,Id=1, Showgrid=false, Showline=false)
////|> Chart.withY_AxisStyle("children N", Side=StyleParam.Side.Right,Id=2, Overlaying=StyleParam.AxisAnchorId.Y 1, Showgrid=false, Showline=true)
////|> Chart.withX_AxisStyle("pre_k", Showgrid=false, Showline=true)
////|> Chart.withLegend false
////|> Chart.Show

////let recordPoints = [1.;24.;25.;26.;28.;32.]

////let node =
////    best_KKZ_det |> Auxilliary.Tree.findNode ["35";"2";"|p104|p111|p112|p16|p25|p32|p53|p72|p88|p103|p110|p31|p36|p44|p78|p89|p118|p119|p19|p21|p56|p7|p90|p94|p116|p3|p42|p11|p70|p86|p45|p67|p114|p13|p33|p77|p105|p117|p18|p38|p41|p55|p59|p6|p65|p71|p85|p91|p99|p23|p60"]
////    //best_KKZ_det |> Auxilliary.Tree.findNode ["35";"2";"p115|p2|p46|p5|p62|p69|p86|p87"]
////    //best_KKZ_det |> Auxilliary.Tree.findNode ["35";"2";"|p107|p75|p82|p100|p108|p15|p28|p37|p48|p51|p57|p61|p63|p64|p83|p95|p0|p39|p74|p96|p109|p47|p50|p68|p102|p120|p17|p22|p58|p66|p8|p81|p84|p1|p49|p97|p14|p26|p29|p34|p43"]
////    //best_KKZ_det |> Auxilliary.Tree.findNode ["35";"2";"|p23|p60|p70|p104|p111|p112|p16|p25|p32|p53|p72|p88|p103|p110|p31|p36|p44|p78|p89|p105|p117|p18|p38|p41|p55|p59|p6|p65|p71|p85|p91|p99|p118|p119|p19|p21|p45|p67|p90|p94|p114|p13|p33|p77|p11|p116|p3|p42|p56|p7"]
////    //|> Array.concat

////let leaves = (gvsprek.[27-2] |> snd) |> Tree.filterLeaves
////leaves.[4].Length
////Plots.drawKinetik leaves.[4] (recordPoints |> List.toArray)  "path 35 other 2" |> Chart.Show

////let drawLeaves title (tree: Types.Node<string,Types.Item>) =
////    tree
////    |> Auxilliary.Tree.filterLeaves
////    |> (fun x -> x.[0 .. ])
////    |> Array.mapi (fun c list ->
////        let col =
////            match c with
////            |0 -> "rgba(0,0,0,1)"  
////            |1 -> colorOrange
////            |2 -> colorBlue   
////            |3 -> colorGreen 
////            |4 -> colorYellow 
////            |5 -> "rgba(0,180,0,1)"
////            |6 -> "rgba(0,180,180,1)"
////            |7 -> "rgba(180,0,180,1)"
////            |8 -> "rgba(180,180,0,1)"
////            |9 -> "rgba(0,0,180,1)"
////            |10 -> "rgba(180,0,0,1)"
////            |_ -> "rgba(0,0,0,1)"
////        list |> Array.map (fun d -> Chart.Line(recordPoints, d.dataL, sprintf "%i" d.ID, Color=col) ) |> Chart.Combine
////                )
////    |> Chart.Combine
////    |> Chart.withTitle title
////    |> Chart.withSize (600.,400.)
////    |> Chart.Show

////drawLeaves "" (gvsprek.[26-2] |> snd)

/////// draw proteins kinetic as a rangePlot with mean line, for the set except given special lines (optional), 

////drawKinetikRangeStack (recordPoints |> Array.ofList) "" (gvsprek.[27-2] |> snd |> Tree.filterLeaves) |> Chart.withLegend (false) |> Chart.Show


/////
//let data1 = 
//    ChlamyProteome.dataAll
//    |> Array.filter (fun x -> x.BinL.[0]="1")
//    |> Array.mapi (fun id x -> {x with ID=id})

///// SSN with pure combinatoric approach (ideal case, the biggest Gain)
//let treeCombi = Functions.SSN.createTree (f 120) (None) Types.Mode.SSN_combi data1

///// SSN with KMeanSwap with deterministic centroid initialization (KKZ) and without any pre-clustering
//let treeKMS_KKZ = Functions.SSN.createTree (f 120) (None) (Types.Mode.SSN (kMeanSwapDet 120)) data1
//let treeKMS_KKZ_ = Functions.SSN.createTree (f 120) (None) (Types.Mode.SSN (kMeanSwapDetCentr 120)) data1

///// compare two networks

//GePhi.sendToGephiFromTreeParam treeCombi
//GePhi.sendToGephiFromTreeParam treeKMS_KKZ
//GePhi.sendToGephiFromTreeParam treeKMS_KKZ_

//treeCombi.GroupGain  // 2.0381 for 1.3 (3 children)
//treeKMS_KKZ.GroupGain // 1.8467 for 1.3 (4 children)
//treeKMS_KKZ_.GroupGain // 1.9497 for 1.3 (4 children)

//let comparison1 = Tree.treeComparison treeCombi treeKMS_KKZ
//let comparison2 = Tree.treeComparison treeCombi treeKMS_KKZ_

//Plots.drawKinetik ([|data1 |> Array.find (fun x -> x.ID = 14);data1 |> Array.find (fun x -> x.ID = 44)|]) [|1 .. 6|]  "path 1.3" |> Chart.Show

////////
//let data =
//    [|
//    [|0.;0.|];
//    [|1.;0.|];
//    [|0.;1.|];
//    [|1.;1.|];
//    [|2.;2.|];
//    [|3.;2.|];
//    [|2.;3.|];
//    [|3.;3.|];
//    [|2.5;2.5|];
//    [|2.;2.5|];
//    |]

//let matrix = SSN.distanceMatrixWeighted None data
//let nR = 10.
//let labeledDataTrue = Array.zip [|1;1;1;1;2;2;2;2;2;2|] [|0 .. 9|] 
//let labeledDataWrong = Array.zip [|1;1;1;1;1;1;2;2;2;2|] [|0 .. 9|] 

//let parents = [|0 .. 9|]
//let label1 = [|0 .. 3|]
//let label2 = [|4 .. 9|]

///// G=dD*IC, the bigger the better
//let dDistIC (labeledID: (int*int) []) =
//    let labels = labeledID |> Array.map fst |> Array.groupBy (fun i -> i) |> Array.map (fun (l,a) -> (l,a.Length))
//    labeledID
//    |> Array.sumBy (fun (label,id) -> 
//                    let nC = labels |> Array.find (fun (name,_) -> name=label) |> snd |> float
//                    let group = if (label=1) then label1 else label2
//                    let ic = -((nC/nR)*log2(nC/nR)+((nR-nC)/nR)*log2((nR-nC)/nR))
//                    let dDist = (SSN.dSumFn [|id|] parents matrix) - (SSN.dSumFn [|id|] group matrix)
//                    dDist*ic)

///// G= D*IC, the bigger the worse
//let distIC (labeledID: (int*int) []) =
//    let labels = labeledID |> Array.map fst |> Array.groupBy (fun i -> i) |> Array.map (fun (l,a) -> (l,a.Length))
//    labeledID
//    |> Array.sumBy (fun (label,id) -> 
//                    let nC = labels |> Array.find (fun (name,_) -> name=label) |> snd |> float
//                    let group = if (label=1) then label1 else label2
//                    let ic = -((nC/nR)*log2(nC/nR)+((nR-nC)/nR)*log2((nR-nC)/nR))
//                    let dDist = (SSN.dSumFn [|id|] group matrix)
//                    (dDist)/ic)

///// G=dD, the bigger the better
//let dDist (labeledID: (int*int) []) =
//    labeledID
//    |> Array.sumBy (fun (label,id) -> 
//                    let group = if (label=1) then label1 else label2
//                    let dDist = (SSN.dSumFn [|id|] parents matrix) - (SSN.dSumFn [|id|] group matrix)
//                    dDist)

///// G= D, the bigger the worse
//let dist (labeledID: (int*int) []) =
//    labeledID
//    |> Array.sumBy (fun (label,id) -> 
//                    let group = if (label=1) then label1 else label2
//                    let dDist = (SSN.dSumFn [|id|] group matrix)
//                    dDist)

///// testing kMS

//let testDDICt = dDistIC labeledDataTrue
//let testDDICf = dDistIC labeledDataWrong


//let testDICt = distIC labeledDataTrue
//let testDICf = distIC labeledDataWrong


//let testDDt = dDist labeledDataTrue
//let testDDf = dDist labeledDataWrong


//let testDt = dist labeledDataTrue
//let testDf = dist labeledDataWrong

//// ok, and now I have to get numbers from theoretical calculations...