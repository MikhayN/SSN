// Learn more about F# at http://fsharp.org
// See the 'F# Tutorial' project for more help.

open FSharp.Plotly
open FSharpAux
open System 
open System.IO

[<EntryPoint>]
let main argv =
    printfn "start %A" argv

    let path = "1"

    printfn "initialized path %s" path

    let data =
        TestData.ChlamyProteome.dataAll
        |> Array.filter (fun x -> x.BinL.[0]=path)
        |> Array.mapi (fun id x -> {x with ID=id})

    printfn "read data, size = %i" data.Length

    let x = Functions.applySST_walkFromHC (data.Length) data

    printfn "tree is created, GroupGain = %f" x.GroupGain

    //let items1 = Auxilliary.Tree.findNodeMembers ["0";"mix-p0-p33-p6-p11-p36-p4-p14-p15-p16-p8-p31-p19-p3-p30-p21-p22-p35-p34-p37-p26"] x
    //let items2 = Auxilliary.Tree.findNodeMembers ["0";"mix-p10-p1-p29-p38-p40-p32-p39-p28"] x
    //let items3 = Auxilliary.Tree.findNodeMembers ["0";"mix-p13-p23-p18-p24-p2-p5-p25"] x

    //Plots.drawKinetik items1 [|1. .. 1. .. 6.|] "" |> Chart.Show
    //Plots.drawKinetik items2 [|1. .. 1. .. 6.|] "" |> Chart.Show
    //Plots.drawKinetik items3 [|1. .. 1. .. 6.|] "" |> Chart.Show

    //Auxilliary.Write.allPathsWrite [|"0"|] TestData.ChlamyProteome.dataAll |> ignore

    let treeLinesW =
        x
        |> Auxilliary.Tree.filterLeaves
        |> Array.concat
        |> Array.map (fun x ->  (sprintf "%i" x.ID) + "\t" + (String.Join(";", x.OriginalBin))  + "\t" + (String.Join(";", x.BinL)))

    printfn "tree is flattened, size = %i" treeLinesW.Length

    File.AppendAllLines(sprintf @"c:\Users\mikha\source\repos\SSN\results\chlamy_path%s_debugPurpose.txt" path, treeLinesW)

    0 // return an integer exit code
