module Plots

open Types
open Auxilliary

open System
open FSharp.Plotly
open FSharpAux
open FSharp.Stats

/// colours
let colorBlue = "rgba(91,155,213,1)"
let colorOrange = "rgba(237,125,49,1)"
let colorGray = "rgba(165,165,165,1)"
let colorYellow = "rgba(255,192,0,1)"
let colorGreen = "rgba(112,173,71,1)"
let colorGrayBlue = "rgba(68,84,106,1)"
let colorBrightBlue = "rgba(68,114,196,1)"

/// draw a histogram, with given bin (optional) and a vertical line, defined by position on 0X and label (optional)
let drawHistogram name (bin: float option) (point: (float*string) option) data =
    let bw = 
        match bin with
        |Some x -> x
        |None -> FSharp.Stats.Distributions.Bandwidth.nrd0 (data |> Array.ofList)

    let temporal = 
        data
        |> Distributions.Frequency.create bw 
        |> Distributions.Frequency.getZip

    let drawPoint =
        match point with
        |Some (v,anno) -> 
                let max =  
                    temporal
                    |> Seq.toList
                    |> List.maxBy snd 
                ([(v, float (snd max));(v, 0.0)], anno)
        |None -> ([],"")

    let hist = Chart.Column (temporal, name)
    let line = Chart.Line (fst drawPoint, snd drawPoint, Labels = [snd drawPoint; ""])
    [hist;line] 
    |> Chart.Combine 
    |> Chart.withTraceName name


/// draw proteins kinetic as lines
let drawKinetik (data: Item []) time title =
    
    let dataLine anno (protein: Item)  =
        let data = protein.dataL
        let points = Array.zip time data
        Chart.Line (points, Name = (sprintf "%i" protein.ID), Labels = [anno;"";"";"";"";""])
        
    data 
    |> Array.map (fun i -> dataLine "" i)
    |> Chart.Combine
    |> Chart.withTitle title
//    |> Chart.Show

/// draw proteins kinetic as lines
let drawKinetikTitle (data': Node<string,Item>) (binList : string list) time =
    let data = data' |> Tree.findNodeMembers binList
    let dataLine anno (protein: Item)  =
        let data = protein.dataL
        let points = Array.zip time data
        Chart.Line (points, Name = (sprintf "%i" protein.ID), Labels = [anno;"";"";"";"";""])
         
    data 
    |> Array.map (fun i -> dataLine "" i)
    |> Chart.Combine
    |> Chart.withTitle (String.Join(".", binList))
//    |> Chart.withSize (600,400)
//    |> Chart.ShowWPF

/// draw proteins kinetic as lines with given length
let drawKinetikTime (data: Item []) (timeLine: float [])=
    
    let dataLine anno (protein: Item)  =
        let data = protein.dataL
        let points = Array.zip timeLine data
        Chart.Line (points, Name = (sprintf "p%i: bin%A" protein.ID protein.BinL), Labels = [anno;"";"";"";"";""])
    
    data 
    |> Array.map (fun i -> dataLine "" i)
    |> Chart.Combine


/// draw proteins kinetic as a rangePlot with mean line, for the set except given special lines (optional), 
///specified by index in a set of items and a label for the line
let drawKinetikRange time title (data: Item [] [])  =

    let dataLine label (protein: Item)  =
        let data = protein.dataL
        let points = Array.zip time data
        let name = string protein.ID
        Chart.Line (points, Name = name)//, Labels = [label;"";"";"";"";""])
    
    let dataLineAsRange anno (proteins: Item [])  =
        let (mean,min,max) = 
            proteins 
            |> Array.map (fun protein -> protein.dataL) 
            |> JaggedArray.transpose 
            |> Array.map (fun timepoint -> (Array.average timepoint, Array.min timepoint, Array.max timepoint))
            |> Array.unzip3

        let meanTime = Array.zip time mean 
        let rangePlot = Chart.Range (meanTime, min, max)//, Labels = [anno;"";"";"";"";""])
        rangePlot //meanLine

    let specialData = 
        data
        |> Array.filter (fun i -> i.Length=1)
        |> Array.concat

    let range =
        data
        |> Array.filter (fun i -> i.Length>1)
        |> Array.map (dataLineAsRange "")

    let specialLines =
        specialData
        |> Array.map (fun protein -> dataLine (string protein.ID) protein)

    [range; specialLines] 
    |> Array.concat 
    |> Chart.Combine 
    |> Chart.withTraceName title 

///specified by index in a set of items and a label for the line
let drawKinetikRangeStack time title (data: Item [] [])  =

    let minY = data |> Array.map (Array.map (fun i -> i.dataL |> Array.min) >> Array.min) |> Array.min
    let maxY = data |> Array.map (Array.map (fun i -> i.dataL |> Array.max) >> Array.max) |> Array.max

    let dataLine col (protein: Item)  =
        let data = protein.dataL
        let points = Array.zip time data
        let name = string protein.ID
        Chart.Line (points, Name = name, Color=col)
    
    let dataLineAsRange col (proteins: Item []) =
        let (mean,min,max) = 
            proteins 
            |> Array.map (fun protein -> protein.dataL) 
            |> JaggedArray.transpose 
            |> Array.map (fun timepoint -> (Array.average timepoint, Array.min timepoint, Array.max timepoint))
            |> Array.unzip3

        let meanTime = Array.zip time mean 
        let rangePlot = 
            Chart.Range (time, mean, min, max, Name = sprintf "%s-%A" (string proteins.[0].ID) (proteins.[0].BinL.[proteins.[0].BinL.Length-1]) , Color=col) 
            |> Chart.withLineStyle (Color=col, Dash=StyleParam.DrawingStyle.Solid)
        let meanLine = Chart.Line (time, mean, Color=col)
        [rangePlot; meanLine] 
        |> Chart.Combine 
        |> Chart.withY_AxisStyle ("", (minY,maxY), Showline=true, Showgrid=false) 
        |> Chart.withX_AxisStyle ("", (time.[0],time.[time.Length-1]), Showline=false, Showgrid=false)

    let singletons =
        data
        |> Array.filter (fun i -> i.Length=1)
        |> Array.mapi (fun c protein -> 
            let col =
                match c with
                |0 -> "rgba(0,0,0,1)"   
                |1 -> colorOrange
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
            dataLine col protein.[0]
        )
        |> Chart.Combine 
        |> Chart.withY_AxisStyle ("", (minY,maxY), Showline=true, Showgrid=false)
        |> Chart.withX_AxisStyle ("", (time.[0],time.[time.Length-1]), Showline=false, Showgrid=false)

    let range =
        data
        |> Array.filter (fun i -> i.Length>1)
        |> Array.mapi (fun c list ->
            let col =
                match c with
                |0 -> "rgba(0,0,0,1)"
                |1 -> colorOrange  
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
            dataLineAsRange col list
            )

    [|range;[|singletons|]|] 
    |> Array.concat 
    |> Chart.Stack 2
    //|> Chart.Combine 
    |> Chart.withTraceName title 

let drawLeaves title recordPoints (tree: Types.Node<string,Types.Item>) =
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
