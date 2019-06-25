module Types

/// basic type for data (after processing)
type Item = {
    ID          : int
    ProteinL    : string array
    OriginalBin : string array
    BinL        : string array 
    dataL       : float array
    }

/// Generic tree node containing member list and child map
type Node<'key, 'data when 'key : comparison> = 
    {
    GroupGain : float     //group gain for the node: max from stepGain and confGain
    StepGain : float      //step gain for the node 
    ConfGain : float*(string list)       //conf gain for the node
    Member   : 'data array
    Children : Map<'key, Node<'key,'data>> 
    }
    
type KMeanSwap = (float -> float -> int -> int -> float) -> float [,] -> int -> Map<string, Item []> -> Map<string, Item []> []
type PreCluster = int -> seq<float> option -> Item list -> Map<string, Item []>


type ClusterFn = int -> Map<string, Item []> -> (string * (Item [])) [] []
type WalkingFn = int -> float [,] -> Map<string,Node<string,Item>> -> (float -> float -> int -> int -> float) -> Map<string, Item []> -> Map<string, Item []>
type OnlyClusterFn = int -> float [,] -> Map<string,Node<string,Item>> -> (float -> float -> int -> int -> float) -> Map<string, Item []> -> Map<string, Item []>

/// creating/alterating tree approach
type Mode =   
    |MM_raw                                 // MapMan only annotated subbins
    |MM                                     // MapMan with broken leaves
    |SSN_clustering of OnlyClusterFn
    |SSN of (KMeanSwap)                     // with KMean Swap (scheme-wise approximation)
    |SSN_pre of (KMeanSwap * PreCluster)    // with KMean Swap (scheme-wise approximation) and pre-clustering (k_max decreasing)
    |SST_walk of (WalkingFn)    // with kMean as a start point for gain walking
    |SSN_combi                              // without simplification, pure combinatorics


