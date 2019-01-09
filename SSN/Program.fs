// Learn more about F# at http://fsharp.org
// See the 'F# Tutorial' project for more help.

open FSharp.Plotly

[<EntryPoint>]
let main argv =
    printfn "start %A" argv
    let x = FunctionsExp.applySSNold TestData.ChlamyProteome.data 100
    GePhi.sendToGephiFromTreeParam x |> ignore

    let items1 = Auxilliary.Tree.findNode ["0";"mix-p0-p33-p6-p11-p36-p4-p14-p15-p16-p8-p31-p19-p3-p30-p21-p22-p35-p34-p37-p26"] x
    let items2 = Auxilliary.Tree.findNode ["0";"mix-p10-p1-p29-p38-p40-p32-p39-p28"] x
    let items3 = Auxilliary.Tree.findNode ["0";"mix-p13-p23-p18-p24-p2-p5-p25"] x

    Plots.drawKinetik items1 [|1. .. 1. .. 6.|] "" |> Chart.Show
    Plots.drawKinetik items2 [|1. .. 1. .. 6.|] "" |> Chart.Show
    Plots.drawKinetik items3 [|1. .. 1. .. 6.|] "" |> Chart.Show

    Auxilliary.Write.allPathsWrite [|"0"|] TestData.ChlamyProteome.data |> ignore

    0 // return an integer exit code
