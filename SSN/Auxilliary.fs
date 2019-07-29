module Auxilliary

open Functions
open Functions.General

open Types
open TestData.General

open System
open System.IO

open FSharp.Plotly
open FSharpAux

module Tree =

    /// Generic tree node empty
    let emptyP = { Member = [||] ; GroupGain = 0.; StepGain = 0.; ConfGain = 0.,[]; Children = Map.empty }

    ///
    let mapTreeToList f tree' =
        let rec loop depth tree =
            [match tree.Children with
                | x when x = Map.empty -> 
                        yield (f tree)
                | _  -> 
                        yield (f tree)
                        for child in tree.Children do                                        
                            yield! loop (depth+1) child.Value 
            ]
        loop 0 tree'

    /// take a tree and return a list of (depth*(node as a tree)) -> very memory-consuming!!!
    let mapTreeFlat tree' =
        let rec loop depth tree =
            [|match tree.Children with
                | x when x = Map.empty -> 
                        yield (depth, tree)
                | _  -> 
                        yield (depth, tree)
                        for child in tree.Children do                                        
                            yield! loop (depth+1) child.Value 
            |]
        loop 0 tree'

    /// take a tree and flat it by return a list of (depth*(node with empty children tree)
    let mapTreeFlat' tree' =
        let rec loop depth tree =
            [match tree.Children with
                | x when x = Map.empty -> 
                        yield (depth, tree)
                | _  -> 
                        yield (depth, {tree with Children=(Map.ofList ["",emptyP])})
                        for child in tree.Children do                                        
                            yield! loop (depth+1) child.Value 
            ]
        loop 0 tree'

    /// take a tree and return a list of (depth*(members of node))
    let mapTreeMembers tree' =
        let rec loop depth tree =
            [match tree.Children with
                | x when x = Map.empty -> 
                        yield (depth, tree.Member)
                | _  -> 
                        yield (depth, tree.Member)
                        for child in tree.Children do                                        
                            yield! loop (depth+1) child.Value 
            ]
        loop 0 tree'

    /// map node bins
    let mapTreeBins tree' =
        let rec loop depth tree acc =
            [match tree.Children with
                | x when x = Map.empty -> 
                        yield (depth, acc)
                | _  -> 
                        yield (depth, acc)
                        for child in tree.Children do                                        
                            yield! loop (depth+1) child.Value (child.Key::acc)
            ]
        loop 0 tree' []

    let mapTreeToBinsMembers tree' =
        let rec loop tree acc =
            [|match tree.Children with
                | x when x = Map.empty -> 
                        yield (acc |> List.rev, tree.Member)
                | _  -> 
                        //yield (depth, acc, tree.Member)
                        for child in tree.Children do                                        
                            yield! loop child.Value (child.Key::acc)
            |]
        loop tree' []
        |> Array.map (fun (bin,xList) -> xList |> Array.map (fun x -> bin,x) )
        |> Array.concat

    /// find node as a tree with given keyPath 
    let findNode (keyPath: 'a list) tree = 
        let rec loop (keyPath: 'a list) tree = 
            match keyPath with
            | [] -> tree
            | k::key ->             
                match Map.tryFind k tree.Children with
                | Some tree -> loop key tree
                | None -> emptyP
        loop keyPath.Tail tree

    /// find nodeMembers with given keyPath 
    let findNodeMembers (keyPath: 'a list) tree = 
        let rec loop (keyPath: 'a list) tree = 
            match keyPath with
            | [] -> tree.Member
            | k::key ->             
                match Map.tryFind k tree.Children with
                | Some tree -> loop key tree
                | None -> [||]
        loop keyPath.Tail tree

    /// find a node with given size
    let rec filterSizesList tree =
        [match tree.Children with
                | x when x = Map.empty -> 
                        yield tree.Member.Length
                | _  -> 
                        yield tree.Children.Count
                        for child in tree.Children do                                        
                            yield! filterSizesList child.Value 
        ]

    /// count number of children for each node in a tree (0's are not shown)
    let rec filterChildrenList tree =
        [match tree.Children with
                | x when x = Map.empty -> 
                        yield 0
                | _  -> 
                        yield tree.Children.Count
                        for child in tree.Children do                                        
                            yield! filterChildrenList child.Value 
        ]
        |> List.filter (fun i -> i <> 0)

    /// filter all leaves as flatted tree
    let rec filterLeaves tree =
        [|match tree.Children with
                | x when x = Map.empty -> 
                        yield tree.Member
                | _  -> 
                        for child in tree.Children do                                        
                            yield! filterLeaves child.Value 
        |]

    /// count pairs, not in the same clusters between trees. It reqiures all items to be identical in both trees 
    let treeComparison (treeA: Node<string,Item>) (treeB: Node<string, Item>) : int =
        let m tree = // adjucency matrix for the tree
            let idBin = 
                tree 
                |> filterLeaves |> Array.concat |> Array.map (fun x -> x.ID,x.BinL) |> Array.sortBy fst
            idBin |> Array.allPairs idBin |> Array.map (fun ((id1,bin1),(id2,bin2)) -> if bin1=bin2 then 1 else 0)
        let lA = m treeA 
        let lB = m treeB 
        Array.zip lA lB |> Array.sumBy (fun ((a),(b)) -> if a=b then 0 else 1)


module ClusterCheck =

    let checkSgain f (matrix: float [,]) (conf: Item [] []) =
        let parentGroup = conf |> Array.concat
        let gFn (current: Item []) = getStepGainFn f (groupIDFn current) (groupIDFn parentGroup) parentGroup.Length matrix
        conf |> Array.map (gFn) |> Array.sum

    let checkSilhouette (matrix: float [,]) (config: Item [] []) =
        
        let silhouette (item:int) (conf: int [] []) (matrix: float [,]) =
            let (eigenCluster,others) = conf |> Array.partition (fun cls -> cls |> Array.contains item)
            let ai = eigenCluster.[0] |> Array.averageBy (fun j -> matrix.[item,j])
            let bi = if others = [||] then ai else others |> Array.map (fun cl -> cl |> Array.averageBy (fun j -> matrix.[item,j])) |> Array.min
            (bi-ai)/(max ai bi)

        let idList = config |> Array.concat |> Array.map (fun p -> p.ID)
        let conf = config |> Array.map (fun list -> list |> groupIDFn)
        idList |> Array.averageBy (fun i -> silhouette i conf matrix)

module Analysis =

    let findMaxDistInMatrix idCurrent (idGroup: int array) (matrix: float [,]) = 
            idGroup
            |> Array.fold (fun maxSoFar i -> max maxSoFar matrix.[i,idCurrent]) 0.

    let pointDxC (matrix: float [,]) (nodes: Item [] []) =
        let nTotal = sqrt (float (matrix.Length))
        let nCluster = float nodes.Length
        let complexity = 
            nCluster/nTotal
        let dissimilarityFn node =
            node
            |> Array.map (fun i -> findMaxDistInMatrix i.ID (groupIDFn node) matrix)
            |> Array.max
        let dissimilarity =
            nodes
            |> Array.sumBy (fun i -> dissimilarityFn i)
            |> (fun x -> x/nCluster)
        (dissimilarity,complexity)

    let leavesAtDepth (depth: int) (flatTree: (int*Node<string, Item>) []) =
        flatTree
        |> Array.choose (fun (i,node) -> 
            if i=depth then
                Some (node.Member)
            elif (i<depth && node.Children=Map.empty) then
                Some (node.Member)
            else None )

    let SO_points_Fn matrix tree =
        let flatTree = 
            tree
            |> Tree.mapTreeFlat
        let maxDepth =
            flatTree
            |> Array.maxBy (fun (i,node) -> i)
            |> fst
        [for a in [0 .. maxDepth] do 
            yield 
                flatTree
                |> leavesAtDepth a 
                |> pointDxC matrix]

    let drawComparePlot dataDO (dataSO: (float*float) list) =
    
        let doP = 
            Chart.Point ([dataDO], Name = "DO", Labels = ["DO"])
    
        let soP =
            let anno = [0 .. (dataSO.Length-1)] |> List.map string
            Chart.Line (dataSO, Name = "SO", Labels = anno)
   
        [soP;doP]
        |> Chart.Combine 


    let allPaths (itemsAll: Item []) f weight listN =

        listN
        |> List.map (fun n ->

                        let itemsOfN = 
                            itemsAll 
                            |> Array.filter (fun i -> i.BinL.Length>0 && i.BinL.[0] = (string n))
                            |> Array.mapi (fun index i -> {i with ID=index; dataL = zScoreTransform i.dataL} )

                        let matrixOfN = 
                            itemsOfN
                            |> distMatrixWeightedOf distanceMatrixWeighted (Some weight)
                            /// normalised
                        let dissMax = 
                            matrixOfN
                            |> Array2D.array2D_to_seq
                            |> Seq.max

                        let normMatrix = 
                            matrixOfN
                            |> Array2D.map (fun i -> i/dissMax)

                        let treeBrokenSOofN = readMM  itemsAll.Length  itemsOfN
                        let fastTreeOfN = applySSN  itemsAll.Length itemsOfN

                        let DO_point_N = 
                            fastTreeOfN
                            |> Tree.filterLeaves
                            |> pointDxC normMatrix
    
                        let SO_points_N =
                            SO_points_Fn normMatrix treeBrokenSOofN

                        (string n,DO_point_N::SO_points_N))

    let drawComparePlot3D (data: (string*(float*float) list) list) =
    
        let data3d (s,data') =
            data'
            |> List.map (fun (x,y) -> (s,x,y))

        let onePath s (data': (string*float*float) list) =
            let anno = 
                [for a in [0 .. (data'.Tail.Length-1)] do yield string a]
            [
            Chart.Scatter3d 
                ([data'.Head], StyleParam.Mode.Markers_Text, Name = (sprintf "DO_of_%s" s), 
                Labels = [sprintf "DO_%s" s], Color = "rgba(50,140,140,1)", TextPosition = StyleParam.TextPosition.BottomCenter);
            Chart.Scatter3d (data'.Tail,StyleParam.Mode.Lines_Markers_Text, Name = s, Labels = anno, Color = "rgba(190,140,40,1)")]

        data
        |> List.map (fun (s,xy) -> onePath s (data3d (s,xy)))
        |> List.concat
        |> Chart.Combine 
        |> Chart.withX_AxisStyle("Path N")
        |> Chart.withY_AxisStyle("Dissimilarity")
        |> Chart.withZ_AxisStyle("Complexity")
        |> Chart.withSize(1000.,1000.)
        |> Chart.Show

    let drawSurface (data: (string*(float*float) list) list) =
        let sortedData =
            data
            |> List.sortBy (fun (path,xyList) -> 
                let area =
                    [1 .. (xyList.Length-2)]
                    |> List.map (fun i -> ((snd xyList.[i])+(snd xyList.[i+1]))*((fst xyList.[i+1])-(fst xyList.[i]))/2.)
                    |> List.sum
                area)
            |> List.mapi (fun i (s,(xy)) -> (i,s,(xy)))
        let grid (data': (float*float) list ) =
            [|0. .. 0.03 .. 1.|]
            |> Array.map (fun x0 -> 
                let (x1,y1) = 
                    data'
                    |> List.rev
                    |> List.filter (fun (x,y) -> x<=x0)
                    |> List.last
                let (x2,y2) =
                    data'
                    |> List.rev
                    |> List.filter (fun (x,y) -> x>x0)
                    |> List.head 
                (y1-y2)/(x1-x2)*x0+(y2-(y1-y2)/(x1-x2)*x2)) // something....
        let surface =
            Chart.Surface ((sortedData |> List.map (fun (i,s,xy) -> grid xy.Tail)), Opacity=0.8, Colorscale=StyleParam.Colorscale.Hot)
        let line =
            let points =
                sortedData
                |> List.map (fun (i,s,xy) -> ((fst xy.Head)*30.,i,(snd xy.Head)))
            let labels =
                sortedData
                |> List.map (fun (i,s,xy) -> s)
            Chart.Scatter3d (points,StyleParam.Mode.Lines_Markers,Color="#2ca02c",Labels=labels)
        [surface;line]
        |> Chart.Combine
        |> Chart.withX_AxisStyle("Dissimilarity(*30)")
        |> Chart.withY_AxisStyle("Path")
        |> Chart.withZ_AxisStyle("Complexity")
        |> Chart.withSize(1000.,1000.)
        |> Chart.Show

    // Calculate line intersection
    let calcIntersection (a:(float*float)) (b:(float*float)) (c:(float*float)) (d:(float*float)) =
        let (Ax,Ay),(Bx,By),(Cx,Cy),(Dx,Dy) =
            (a,b,c,d)
        let d = (Bx-Ax)*(Dy-Cy)-(By-Ay)*(Dx-Cx)  

        if  d = 0. then
        // parallel lines ==> no intersection in euclidean plane
            None
        else
            let q = (Ay-Cy)*(Dx-Cx)-(Ax-Cx)*(Dy-Cy) 
            let r = q / d
            let p = (Ay-Cy)*(Bx-Ax)-(Ax-Cx)*(By-Ay)
            let s = p / d

            if r < 0. || r > 1. || s < 0. || s > 1. then
                None // intersection is not within the line segments
            else
                Some((Ax+r*(Bx-Ax)), (Ay+r*(By-Ay)))  // Px*Py

    let intersectDOandSOlines (doPoint: float*float) (soPoints: (float*float) list) =
        let (a),(b) = ((0.,0.),((fst doPoint)/(snd doPoint),1.))
        soPoints
        |> List.pairwise
        |> List.map (fun (x,y) -> calcIntersection x y a b)

    /// return two values: (dist to SO line * dist to DO point)
    let getDistFromZero doPoint =
        (intersectDOandSOlines doPoint)
        >> List.find (fun i -> Option.isSome i)
        >> Option.get
        >> (fun (x,y) -> 
                            (
                            weightedEuclidean None [x;y] [0.;0.],                       // dist to SO line
                            weightedEuclidean None [fst doPoint;snd doPoint] [0.;0.])   // dist to DO point
                            )

    ////
    
    /// draw dissimilarity VS complexity plot
    let drawDCsinglePath weightL items soTree doTree =
        let matrixPath = 
            items
            |> distMatrixWeightedOf distanceMatrixWeighted (weightL)

        let pointDO = 
            pointDxC matrixPath (doTree |> Tree.filterLeaves)

        let pointsSO =
            SO_points_Fn matrixPath soTree

        let normalize (maxDissim: float) (input: float*float) =
            let (x,y) = input
            (x/maxDissim,y)
    
        let maxDiss = 
            pointsSO
            |> List.map fst
            |> List.max

        let pointDOnorm = 
            normalize maxDiss pointDO  

        let pointSOnorm = 
            pointsSO
            |> List.map (fun i -> normalize maxDiss i) 

        let intersectsDOxSO = intersectDOandSOlines pointDOnorm pointSOnorm
    
        [drawComparePlot pointDOnorm (pointSOnorm);
        Chart.Line [(0.0, 0.0); (intersectsDOxSO |> List.choose (fun i -> i) |> List.head)];
        Chart.Line [(0.0, 0.0); pointDOnorm];
        ]
        |> Chart.Combine
        |> Chart.withSize (500.,500.)
        |> Chart.Show

    let getDCmeasure weightL items soTree doTree =
        let matrixPath = 
            items
            |> distMatrixWeightedOf distanceMatrixWeighted (weightL)

        let pointsSO =
            SO_points_Fn matrixPath soTree

        let normalize (maxDissim: float) (input: float*float) =
            let (x,y) = input
            (x/maxDissim,y)
    
        let maxDiss = 
            pointsSO
            |> List.map fst
            |> List.max 

        doTree 
        |> Tree.filterLeaves
        |> pointDxC matrixPath 
        |> normalize maxDiss
        |> (fun (x,y) -> sqrt(x*x+y*y))


module Write =

    let allPathsWrite listN items =
        //let pathFile = @"c:\_n_mikhaylenko\Code_FSharp\Projects\SSN\results\check\"
        //let pathFile = @"D:\Nathan\results\tablesAraArraysOnly\"
        let pathFile = @"..\results\"

        listN
        |> Array.mapi (fun i n ->
                    
                        let itemsOfN =
                            items
                            |> Array.filter (fun i -> i.BinL.Length>0 && i.BinL.[0] = n)
                            |> Array.mapi (fun index i -> {i with ID=index; dataL = i.dataL |> zScoreTransform })

                        let treeMMofN = readMM  items.Length itemsOfN
                        let stopwatch = new System.Diagnostics.Stopwatch()
                        stopwatch.Start()
                        let treeSSNofN = Functions.applySSNold  items.Length itemsOfN             //// change here which SSN to use
                        let time = sprintf "SSN tree calculated in %f s since start" (stopwatch.Elapsed.TotalSeconds)
                        stopwatch.Stop()
                        let dc = Analysis.getDCmeasure None itemsOfN treeMMofN treeSSNofN
                        let fileName = n
                        let title = sprintf "bin: %s, items: %i" n itemsOfN.Length
                        let dcText = sprintf "DxC measure: %f" dc
                        let header = "Bin\tItem"
                        let content = treeSSNofN |> Tree.mapTreeToBinsMembers |> Array.toList |> List.map (fun (bin,x) -> sprintf "%s\t%A" (String.Join(".", bin)) x.ProteinL)
                        //treeArray.[i] <- fastTreeOfN
                        File.WriteAllLines((sprintf "%s%s.txt" pathFile fileName), title :: time :: dcText :: header :: content)
                        )

    ////
    //let writeFileFromResult proteinList (fileName: string) =
    //    let proteinAbb =
    //        (File.ReadAllLines (@"c:\_n_mikhaylenko\BioDB\Data\PhytozomeV11\Creinhardtii\annotation\Creinhardtii_281_v5.5.geneName.txt"))
    //        |> Array.map (fun s -> s |> String.toWords |> Seq.toList |> List.cutAfterN 2 |> fst)
    //    let makeLine protein =
    //        let id = protein.ID
    //        let name =  
    //            try
    //                protein.ProteinL |> Array.map (fun protein -> proteinAbb |> Array.find (fun line -> (line.Head)=protein) |> List.tail) |> List.concat
    //            with
    //                |ex -> [""]
 
    //        sprintf "%i\t%A\t%A" id protein.ProteinL name
    //    let header = "Protein_Group_id\tProteins_Code\tProteins_Names"
    //    let lines =
    //        proteinList
    //        |> List.map (fun i -> makeLine i)
        
    //    File.WriteAllLines (sprintf "c:\\_n_mikhaylenko\\Code_FSharp\\Projects\\DynamicOntology\\results\\tables\\%s.txt" fileName,header::lines)

    //let writeFileSMlist (proteinSList: (string list) list) (proteinMList: (string list) list) (tree: Node<string, Item>) (fileName: string) =
    //    let pathAbb = @"c:\_n_mikhaylenko\BioDB\Data\PhytozomeV11\Creinhardtii\annotation\Creinhardtii_281_v5.5.geneName.txt"
    //    let pathFile = "c:\\_n_mikhaylenko\\Code_FSharp\\Projects\\DynamicOntology\\results\\tables\\"
    //    let proteinAbb =
    //        (File.ReadAllLines pathAbb)
    //        |> Array.map (fun s -> s |> String.toWords |> Seq.toList |> fun x -> x.[0 .. 2])
    //    let makeLine protein =
    //        let id = protein.ID
    //        let name =  
    //            try
    //                protein.ProteinL
    //                |> Array.toList
    //                |> List.map (fun protein ->
    //                    proteinAbb
    //                    |> Array.find (fun line -> (line.Head)=protein)
    //                    |> List.tail)
    //                |> List.concat
    //            with
    //                |ex -> ["notFound"]
    //        sprintf "%i\t%A\t%A" id protein.ProteinL name

    //    let mergeHeader = "Merged proteins from bin "
    //    let splitHeader = "Splitted proteins from bin "
    //    let header = "Protein_Group_id\tProteins_Code\tProteins_Names"
    //    let linesSplitted =
    //        proteinSList
    //        |> List.map (fun x -> x,(Tree.findNode x tree))
    //        |> List.map (fun (bin,list) -> ((sprintf "%s%A" splitHeader bin) |> String.filter (fun c -> c<>'\"'))::header::( list |> Array.toList |> List.map makeLine))
    //        |> List.concat
    //    let linesMerged =
    //        proteinMList
    //        |> List.map (fun x -> x,(Tree.findNode x tree))
    //        |> List.map (fun (bin,list) -> ((sprintf "%s%A" mergeHeader bin) |> String.filter (fun c -> c<>'\"'))::header::( list |> Array.toList |> List.map makeLine))
    //        |> List.concat       
    //    File.WriteAllLines((sprintf "%s%s.txt" pathFile fileName), (List.append linesSplitted linesMerged))


