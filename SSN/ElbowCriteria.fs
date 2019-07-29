module ElbowCriteria

open FSharpAux
open FSharp.Plotly
open FSharp.Stats

open Types
open Functions
open Functions.General
open FunctionsExp
open TestData

module ElbowCrit =
    
    
    let plotElbow critTitle (criterion: Item [] [] -> float) (kXclusterings: (int*(Item [] [])) list) =
        let xy = 
            kXclusterings
            |> List.map (fun (k,sch) -> (k, criterion sch))
        let titel = sprintf "Elbow criterion: %s"  critTitle
        Chart.Line (xy, ShowMarkers=true)
        |> Chart.withTitle titel
        |> Chart.Show

    //let plotElbow' critTitle fn (kXclusterings: (int*(Item [] [])) list) =
    //    let xy = 
    //        fn kXclusterings
    //    let titel = sprintf "Elbow criterion: %s"  critTitle
    //    Chart.Line (xy, ShowMarkers=true)
    //    |> Chart.withTitle titel
    //    |> Chart.Show

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

    //let hartiganIndex (kXclusterings: (int*(Item [] [])) list) =
    //    let n = kXclusterings.Head |> snd |> Array.concat |> Array.length
    //    let wk = 
    //        kXclusterings
    //        |> List.map (fun (k,sch) -> (wIndex sch))
    //    [for k in [1 .. (wk.Length-1)] do 
    //                                        let gamma = float (n-k-1)
    //                                        yield (k,gamma*(wk.[k-1]-wk.[k])/wk.[k])
    //    ]

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

    //let sseStep (kXclusterings: (int*(Item [] [])) list) =
    //    let n = kXclusterings.Head |> snd |> Array.concat |> Array.length
    //    let wk = 
    //        kXclusterings
    //        |> List.map (fun (k,sch) -> (sse sch))
    //    [for k in [1 .. (wk.Length-1)] do
    //                                        yield (k,(wk.[k-1]-wk.[k]))
    //    ]

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

    //let dbIndexStep (kXclusterings: (int*(Item [] [])) list) =
    //    let n = kXclusterings.Head |> snd |> Array.concat |> Array.length
    //    let wk = 
    //        kXclusterings
    //        |> List.map (fun (k,sch) -> (dbIndex sch))
    //    [for k in [1 .. (wk.Length-1)] do
    //                                        yield (k,(wk.[k-1]-wk.[k]))
    //    ]

module GetK =

    let vectorRejection vector projectionLine =
        let pNorm = Vector.norm projectionLine
        let pUnit = Vector.scale (1./pNorm) projectionLine
        let vProj = Vector.scale (Vector.dot vector pUnit) pUnit
        Vector.norm (vector - vProj)

    let lineToElbow fn kXclusterings = 
        let xy = kXclusterings |> List.map (fun (k,sch) -> ( float k, fn sch))
        let xmax = xy |> List.map (fst) |> List.max
        let xmin = xy |> List.map (fst) |> List.min
        let ymax = xy |> List.map (snd) |> List.max
        let ymin = xy |> List.map (snd) |> List.min
        xy
        |> List.map (fun (x,y) ->( x,  [(x-xmin)/(xmax-xmin);(y-ymin-ymax)/(ymax-ymin)] |> MatrixTopLevelOperators.vector))

    let optimumFound criterion (data: Item []) = 
        
        let pre_k_max =
            2.**(float data.[0].dataL.Length)
            |> int

        let clusteringSet_fn data =
            let basic = 
                data
                |> Array.map (fun protein -> protein.dataL)
                |> ML.Unsupervised.HierarchicalClustering.generate (weightedEuclidean None) (ML.Unsupervised.HierarchicalClustering.Linker.completeLwLinker)

            [1 .. (min data.Length pre_k_max)] |> List.map (fun k -> 
                (k,basic
                |> ML.Unsupervised.HierarchicalClustering.cutHClust k
                |> List.map (List.map (fun i -> data.[ML.Unsupervised.HierarchicalClustering.getClusterId i]) >> List.toArray)
                |> List.toArray ))

        let clusteringSet = clusteringSet_fn data
        let line = lineToElbow criterion clusteringSet // sse wIndex
        let xy = [1.;-1.] |> MatrixTopLevelOperators.vector
        let distances = 
            line 
            |> List.map (fun (k,i) -> (k,(i.[0], vectorRejection i xy)))
        let k = distances |> List.maxBy (snd >> snd) |> fst |> int
        clusteringSet
        |> List.find (fun (ki,x) -> (ki=k) )