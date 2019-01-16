module SOM

open Functions
open Functions.SSN
open Types

open FSharpAux

type SOMclusters =
    {WeightPatterns : float [] [];
    ResultingGroups : (int*int*(float [])) []}

//    module Functions =
// sigma = 8. ; tau = 100.
let neighFn i winnerIndex currentIndex = // neighbourhood function
    let sigma = 8.0
    let tau = 100.0
    exp( -( float ((winnerIndex - currentIndex)*(winnerIndex - currentIndex)) )/(sigma*exp( -(float i)/tau )))

// etaZero = 0.5 ; tau = 1000.
let rate i = // learning rate
    let etaZero = 0.5
    let tau = 1000.0
    etaZero*exp(-(float i)/tau)
 

/// the main function, kinetic data in dataX are rows
let som weight (dataX: float [,]) k power =
    let d = dataX.[0,0..].Length // dimention of data vector
    let n = dataX.[0 ..,0].Length // amount of elements in data
    let r = System.Random()

    // create weigthing functions randomly
    let weightPatterns = (Array2D.init k d (fun row col -> (1.)/(float (r.Next(1,10)))))

    let findWinner x =
        let distArray = Array.init k (fun i -> weightedEuclidean weight dataX.[x,0..] weightPatterns.[i,0..])
        let min = distArray |> Array.min
        let index = distArray |> Array.findIndex (fun i -> i = (distArray |> Array.min))
        (index, min)

    let updateWeights x (winner: int) i =
        for a in [0 .. k-1] do
            (weightPatterns.[a,0..] <- Array.init d (fun elem -> weightPatterns.[a,elem] + (rate i)*(neighFn i winner a)*(dataX.[x,elem]-weightPatterns.[a,elem])))

    let repeatUpdate nTimes (dataX: float [,]) =
        for time in [1 .. nTimes] do
            let x = r.Next(0,(n-1))
            let winner = findWinner x
            updateWeights x (fst winner) time

    repeatUpdate power dataX
    {WeightPatterns = weightPatterns |> Array2D.toJaggedArray;
    ResultingGroups =
        [|0 .. n-1|]
        |> Array.map (fun proteinN -> (proteinN, fst (findWinner proteinN)))
        |> Array.map (fun (itemi,cl) -> (cl,itemi, dataX |> Array2D.toJaggedArray |> Array.item itemi))}

//let useSOM weight (items: Item []) (power:int) (nGroups: int ) : Map<string,Item []> =
//    let data =
//        items
//        |> Array.map (fun i -> i.dataL)
//        |> Array2D.ofJaggedArray

//    let result = som weight data nGroups (power*data.Length)

//    result.ResultingGroups
//    |> Array.groupBy (fun (cl,_,_) -> cl)
//    |> Array.map (fun (_,list) ->
//                                let members = list |> Array.map (fun (c,id,x) -> items.[id] )
//                                let newBin = members |> Array.fold (fun acc p -> sprintf "%s-p%i" acc p.ID) "cl"
//                                (newBin, members))
//    |> Map.ofArray