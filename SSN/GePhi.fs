module GePhi

open Functions
open Types

open FSharpGephiStreamer
open FSharpAux
open FSharpAux.IO
open FSharp.Stats
open System

module GePhi_type =

    type NodeParam =
        {
        ProteinIDs: string;
        NodeBin: string;
        LevelDepth: int;
        Size: float;
        Elements: int;
        groupGain: float;
        stepGain: float;
        confGain: float;
        X_Position: float;
        Y_Position: float
        }

/// convert tree structure to nodes and edges and send them to GePhi Streaming
let sendToGephiFromTreeParam (tree: Node<string, Item>) =
    
    let rec toMemberSeq key depth (tree: Node<string,Item>) =
        let fillNodeFrom (items: Item []) gGain sGain cGain : GePhi_type.NodeParam =
            {
            ProteinIDs = String.Join(";", (SSN.groupIDFn items));
            NodeBin = String.Join(".",Array.append items.[0].BinL.[0 .. (depth-1)] [|key|]);
            LevelDepth = depth;
            Size = (float items.Length)/2.;
            Elements = items.Length;
            groupGain = gGain;
            stepGain = sGain;
            confGain = cGain;
            X_Position = 0.;
            Y_Position = 0.;
            }
        [ 
        match tree.Children with
        |x when x=Map.empty -> 
                    yield fillNodeFrom tree.Member tree.GroupGain tree.StepGain (fst tree.ConfGain) 
                        
        |_  ->      
                    yield fillNodeFrom tree.Member tree.GroupGain tree.StepGain (fst tree.ConfGain) 
                    for child in tree.Children do                                        
                        yield! toMemberSeq child.Key (depth+1) child.Value 
        ]

    let key = tree.Member.[0].BinL.[0]
    let flatTree = toMemberSeq key 0 tree

    let yList =
        let plus index acc =
            if (flatTree.[index].LevelDepth>flatTree.[index-1].LevelDepth) then acc
            else (acc + 4. + (float flatTree.[index-1].Elements)/2. + (float flatTree.[index].Elements)/2.)
        let rec loop index acc =
            match index with
            |0 -> acc
            |_ -> loop (index-1) (plus index acc) 
        [for index in [0..(flatTree.Length-1)] do yield loop index 0.0]

    let xList =
        let maxDepth = 
            flatTree
            |> List.maxBy (fun i -> i.LevelDepth)
            |> (fun i -> i.LevelDepth)
        let count level =
            flatTree
            |> List.filter (fun i -> i.LevelDepth=level)
            |> List.maxBy (fun i -> i.Elements)
            |> (fun i -> i.Elements)
            |> float
            
        let rec loop index acc = 
            match index with
            |0 -> acc
            |_ -> loop (index-1) (acc + 10.+(count index)/2.+(count (index-1))/2.)
        [for level in [0 .. maxDepth] do yield loop level 0.0]

    let nodes = 
        flatTree
        |> List.mapi (fun i node -> 
                                    (i, {node with 
                                            X_Position = xList.[node.LevelDepth];
                                            Y_Position = yList.[i]
                                            }))
    let edges = 
        flatTree 
        |> List.tail
        |> List.mapi ( fun i node -> 
                                let source = 
                                    i - (flatTree 
                                    |> List.cutAfterN (i+1)
                                    |> fst
                                    |> List.rev
                                    |> List.findIndex (fun nodeBefore -> nodeBefore.LevelDepth < node.LevelDepth)) 
                                let target = i+1
                                let level = node.LevelDepth
                                let gain = node.groupGain
                                (source,target,(level,gain)))

    let converterNode ((index,param) : int*GePhi_type.NodeParam) =
        [
        Grammar.Attribute.Label param.NodeBin;
        Grammar.Attribute.UserDef ("Level", param.LevelDepth);
        Grammar.Attribute.UserDef ("Elements", param.Elements);
        Grammar.Attribute.UserDef ("Group_Gain", System.Math.Round(param.groupGain,4));
        Grammar.Attribute.UserDef ("Step_Gain", System.Math.Round(param.stepGain,4));
        Grammar.Attribute.UserDef ("Conf_Gain", System.Math.Round(param.confGain,4));
        Grammar.Attribute.UserDef ("Protein_IDs", param.ProteinIDs);
        Grammar.Attribute.Size param.Size; 
        Grammar.Attribute.PositionX param.X_Position;
        Grammar.Attribute.PositionY param.Y_Position
        ]

    let converterEdge (edge : int*int*(int*float)) =
        match edge with
        |(_, _, c) -> [Grammar.Attribute.UserDef ("Level", fst c); Grammar.Attribute.UserDef ("Group_Gain", snd c);]

    let postNodes = 
        nodes 
        |> List.map (fun a -> Streamer.addNode converterNode (fst a) a)

    let postEdges = 
        edges 
        |> List.mapi (fun i (source,target,att) -> Streamer.addEdge converterEdge i (source) (target) (source,target,att))

    [
    postNodes;
    postEdges
    ]