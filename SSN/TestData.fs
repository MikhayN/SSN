module TestData

open Types
open FSharp.Stats
open FSharpAux
open FSharpAux.IO

open FSharp.Plotly

open System
open System.IO
 
module General =

    /// change the path to the local adress of the repo
    let pathToData = @"c:\Users\mikha\source\repos\SSN\"

    // some special custom input info
    let weightProt = seq [0.125;0.78;0.0625;0.09375;0.1875;0.125] 
    let weightTran = seq [0.125;0.78;0.78;0.0625;0.09375;0.1875;0.125] 

    let timeSeriesProt = [|1.; 2.; 3.; 26.; 28.; 32.|] 
    let timeSeriesTran = [|0.5; 1.; 2.; 3.; 26.; 28.; 32.|] 
    let timeSeriesArab = [|1. .. 1. .. 10.|] 

    let rca1 = "Cre04.g229300"
    let rmt1 = "Cre16.g661350"
    let psrp1 = "Cre05.g237450"
    let calvinCycleBin = ["1"; "3"]
    let psIIsubunitBin = ["1";"1";"1";"2"]
    let plastidRP = ["29";"2";"1";"1"]
    let cytosolicRP = ["29";"2";"1";"2"]

    let hasHeader = true
    let separatorTab = '\t'
    let separatorComma = ','

    /// transform raw data in proteins experimental data using given transformation function
    let transformKineticData (transformF: float [] -> float []) (item: Item) =
        {item with dataL = 
                    item.dataL 
                    |> transformF }

    let zScoreTransform (data: float []) = 
        let stats = data |> Seq.ofArray //|> Seq.stats  
        let mean = stats |> Seq.mean //SummaryStats.mean stats
        let sd = stats |> Seq.stDev  //SummaryStats.stDevPopulation stats
        data
        |> Seq.map (fun x -> (x-mean)/sd)
        |> Seq.toArray

module SynData =
    let data1 = // 11111-22222
        [|{ID = 0;
            ProteinL = [|"name-1-0"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.594924196; -0.6799640241; -0.0256084639; 0.4184221454; 0.7876903374; 1.094384202|];};
        {ID = 1;
            ProteinL = [|"name-1-1"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.623806497; -0.6918776997; 0.01799425629; 0.485284777; 0.8499604241; 0.9624447397|];};
        {ID = 2;
            ProteinL = [|"name-1-2"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.660086232; -0.6015048645; -0.01292323899; 0.5076795715; 0.7051654513; 1.061669313|];};
        {ID = 3;
            ProteinL = [|"name-1-3"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.679000203; -0.5919325127; 0.01478577775; 0.4538017553; 0.9110590795; 0.8912861036|];};
        {ID = 4;
            ProteinL = [|"name-1-4"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.601316534; -0.737273413; 0.09090209465; 0.4741557857; 0.6788978321; 1.094634234|];};
        {ID = 5;
            ProteinL = [|"name-2-5"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.620908382; 0.6052804092; 0.1000292681; -0.4897408503; -0.7310015607; -1.105475648|];};
        {ID = 6;
            ProteinL = [|"name-2-6"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.732824906; 0.4496575389; -0.006152173633; -0.358682348; -0.8236386405; -0.9940092832|];};
        {ID = 7;
            ProteinL = [|"name-2-7"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.661185627; 0.5405163318; 0.08330491781; -0.4533593937; -0.7449967746; -1.086650709|];};
        {ID = 8;
            ProteinL = [|"name-2-8"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.571046142; 0.6659879663; 0.1188755915; -0.4858593321; -0.7234369518; -1.146613416|];};
        {ID = 9;
            ProteinL = [|"name-2-9"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.67751769; 0.5858771221; -0.1169415339; -0.4590948617; -0.5317098437; -1.155648573|];}|]

    let data2 = // 11111-22222
        [|{ID = 0;
            ProteinL = [|"name-1-0"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.249006804; -1.252379798; 0.307344253; 0.9641723123; 0.4014543971; 0.8284156403|];};
        {ID = 1;
            ProteinL = [|"name-1-1"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.904745432; -0.03545858216; 0.235293155; 0.3717654776; 1.04635999; 0.2867853912|];};
        {ID = 2;
            ProteinL = [|"name-1-2"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.40088378; -0.4039629479; -0.4222492444; 0.1571901489; 0.5206182222; 1.549287601|];};
        {ID = 3;
            ProteinL = [|"name-1-3"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.64537762; -0.3346709741; -0.01776556405; 0.7748367887; 1.256548116; -0.03357074684|];};
        {ID = 4;
            ProteinL = [|"name-1-4"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.711570784; -0.1134732279; -0.3248863839; 0.7981286614; 0.2279026571; 1.123899077|];};
        {ID = 5;
            ProteinL = [|"name-2-5"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.312453878; 0.6927622691; -1.024884286; -0.5294283502; 0.60074619; -1.051649701|];};
        {ID = 6;
            ProteinL = [|"name-2-6"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.041926113; 0.5805275977; 0.6952254848; 0.07572638555; -0.8619569292; -1.531448652|];};
        {ID = 7;
            ProteinL = [|"name-2-7"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|0.9182899219; 1.299063406; -0.3361523136; -0.08224203814; -1.504164967; -0.294794009|];};
        {ID = 8;
            ProteinL = [|"name-2-8"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.797260723; 0.4243921594; -0.2510013552; -0.2785737992; -0.7522928488; -0.9397848787|];};
        {ID = 9;
            ProteinL = [|"name-2-9"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.45649186; 1.009780815; -0.1196038671; -0.7154160467; -0.8499496687; -0.7813030926|];}|]

    let data3 = // 00-11-222222
        [|{ID = 0;
            ProteinL = [|"name-0-0"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-0.0636620859; -0.06119196923; -0.01466371581; 0.08523164606; 0.08380272863; 0.07048339629|];};
        {ID = 1;
            ProteinL = [|"name-0-1"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-0.0355264886; -0.08083434143; 0.0183349666; -0.03903575135; 0.05908404028; 0.07797757442|];};
        {ID = 2;
            ProteinL = [|"name-1-2"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.178029051; -1.187539947; -0.0563792052; 0.4131349135; 0.9320466381; 1.076766652|];};
        {ID = 3;
            ProteinL = [|"name-1-3"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.859759663; 0.3637071728; 0.3278055376; -0.3354497542; 0.5808568543; 0.9228398528|];};
        {ID = 4;
            ProteinL = [|"name-2-4"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.496163239; 0.6758421078; -0.1958494622; -0.5284020891; -0.03861355314; -1.409140242|];};
        {ID = 5;
            ProteinL = [|"name-2-5"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.549159199; 0.4256318769; 0.3750266048; -0.5547714976; -0.4735941401; -1.321452043|];};
        {ID = 6;
            ProteinL = [|"name-2-6"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.193397759; 0.9218648743; 0.4772811654; -0.4754283842; -1.183154346; -0.9339610681|];};
        {ID = 7;
            ProteinL = [|"name-2-7"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.667690621; 0.4968942602; -0.3091856537; -0.001867291717; -0.6452393139; -1.208292622|];};
        {ID = 8;
            ProteinL = [|"name-2-8"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.470419431; 0.4840225401; 0.457031256; -0.4363468966; -1.343758816; -0.6313675145|];};
        {ID = 9;
            ProteinL = [|"name-2-9"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.808877991; 0.4811514528; -0.3870354405; -0.6965223446; -0.3446561006; -0.8618155582|];}|]

    let data4 = //// it has to be adjusted
        [|{ID = 0;
            ProteinL = [|"name-1-0"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.695517053; 0.1977952906; -0.27546189; 0.5012073488; -0.05323207048; 1.325208374|];};
        {ID = 1;
            ProteinL = [|"name-1-1"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.159536014; -1.13053636; 0.1159842631; 0.6251900104; 0.1525542806; 1.39634382|];};
        {ID = 2;
            ProteinL = [|"name-1-2"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.518924248; -0.9475592846; 0.2818354427; 0.4399212567; 0.871192729; 0.8735341041|];};
        {ID = 3;
            ProteinL = [|"name-1-3"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.918262683; 0.09206142105; 0.1295746959; 0.1813512622; 0.5188014429; 0.9964738607|];};
        {ID = 4;
            ProteinL = [|"name-1-4"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|-1.495405434; -0.1714587399; 0.01373039242; -0.4199918722; 0.5844246814; 1.488700972|];};
        {ID = 5;
            ProteinL = [|"name-2-5"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.170334643; 1.092401338; -0.004731922306; -0.9101257786; -0.08223106729; -1.265647213|];};
        {ID = 6;
            ProteinL = [|"name-2-6"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.752335164; 0.2229842076; 0.02571703711; -0.06511967118; -0.9478241777; -0.9880925599|];};
        {ID = 7;
            ProteinL = [|"name-2-7"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.241509293; 0.9578576745; 0.354962629; -0.8194239235; -0.5219383865; -1.212967286|];};
        {ID = 8;
            ProteinL = [|"name-2-8"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.582923146; 0.8160805522; -0.3156114193; -0.4361005577; -0.5221598957; -1.125131826|];};
        {ID = 9;
            ProteinL = [|"name-2-9"|];
            OriginalBin = [|"root"|];
            BinL = [|"root"|];
            dataL =
            [|1.25210773; 1.086679189; -0.2835872432; -1.133427748; 0.01944862816; -0.9412205564|];}|]

module ChlamyProteome =

    type NameConverter() = 
        inherit SchemaReader.Attribute.ConverterAttribute()
        override this.convertToObj = 
            SchemaReader.Converter.Single (fun (strs : string) -> 
                                                (strs |> String.split ';' |> Array.map (fun i -> String.replace "p" "" i)) |> box)

    type IntConverter() = 
        inherit SchemaReader.Attribute.ConverterAttribute()
        override this.convertToObj = 
            SchemaReader.Converter.Single (fun (strs : string) -> 
                                                (strs |> (String.split '.') |> Seq.map (fun i -> i |> int |> string ) |> Seq.toArray) |> box)  

    type DoubleArrayConverter() = 
        inherit SchemaReader.Attribute.ConverterAttribute()
        override this.convertToObj = 
            SchemaReader.Converter.Collection (fun (strs : seq<string>) -> 
                                                    (strs |> Seq.map (fun s -> FSharpAux.String.tryParseFloatDefault nan s) |> Seq.toArray) |> box)
        
            /// type for raw data (directly after reading)
    type ProteinItemRead = {
        [<SchemaReader.Attribute.FieldAttribute("ProteinGroup")>]                                      [<NameConverter>]           ProteinGroup    : string []
        [<SchemaReader.Attribute.FieldAttribute("Bin")>]                                               [<IntConverter>]            Bin             : string []
        [<SchemaReader.Attribute.FieldAttribute([| "TP0HS";"TP24HS";"TP1R";"TP2R";"TP4R";"TP8R" |])>]  [<DoubleArrayConverter>]    Features        : float []
        }

    /// convert protein info from ProteinItemRead (data in Seq) to Proteins (data in List)  
    let proteinsToListFrom id (item: ProteinItemRead) = 
        {
        ID = id;
        ProteinL = item.ProteinGroup;
        OriginalBin = item.Bin;
        BinL = item.Bin;
        dataL = item.Features
        } 

    //// Read Proteins database

    let csvPath = sprintf @"%sdata\HeatShockRecovery_sig_WardLinkage.txt" General.pathToData

    let reader    = new SchemaReader.Csv.CsvReader<ProteinItemRead>(schemaMode=SchemaReader.Csv.SchemaMode.Fill)

    /// Variable, contains all raw data from protein DataBase file in csvPath
    let dataProtein = 
        reader.ReadFile (csvPath, General.separatorTab, General.hasHeader) 
        |> Seq.distinctBy (fun x -> x.ProteinGroup)
        |> Seq.toArray

    //let data = 
    //    dataProtein 
    //    |> Array.filter (fun i -> i.Bin.Length>0 && i.Bin.[0] = "1")
    //    |> Array.mapi (fun i p -> 
    //        proteinsToListFrom i p 
    //        |> General.transformKineticData General.zScoreTransform 
    //        |> fun x -> {x with BinL=[|x.BinL.[0]|]; OriginalBin=[|x.OriginalBin.[0]|]})
    //    |> Array.chunkBySize 50
    //    |> Array.item 0

    let dataAll = 
        dataProtein
        |> Array.mapi (fun i p -> 
            proteinsToListFrom i p 
            |> General.transformKineticData General.zScoreTransform )

module ChlamyTranscriptome =

    type BinConverter() = 
        inherit SchemaReader.Attribute.ConverterAttribute()
        override this.convertToObj = 
            SchemaReader.Converter.Single (fun (strs : string) -> 
                                                (strs |> String.filter (fun c -> c<>''')|> String.split '.' ) |> box)

    type IdentifierConverter() = 
        inherit SchemaReader.Attribute.ConverterAttribute()
        override this.convertToObj = 
            SchemaReader.Converter.Single (fun (strs : string) -> 
                                                (strs |> String.split '|' |> Array.item 0 |> String.replace "p" "") |> box)

    // read MapMan to extract Bin

    type MapManRead = {
        [<SchemaReader.Attribute.FieldAttribute("BINCODE")>]    [<BinConverter>]       BinString           : string []
        [<SchemaReader.Attribute.FieldAttribute("IDENTIFIER")>] [<IdentifierConverter>]       ProteinIdentifier   : string
        }

    let csvPathMM = sprintf @"%sdata\Creinhardtii_236 - ManMapList.txt" General.pathToData

    let readerMM    = new SchemaReader.Csv.CsvReader<MapManRead>(schemaMode=SchemaReader.Csv.SchemaMode.Fill)

    //// Read Proteins database

    type ProteinItemReadT = {
        [<SchemaReader.Attribute.FieldAttribute("Key")>]                                      [<ChlamyProteome.NameConverter>]           ProteinGroup    : string []
        [<SchemaReader.Attribute.FieldAttribute(
            [| "HS.1-HS.0.logFC";"HS.2-HS.0.logFC";"HS.3-HS.0.logFC";
            "R.1-HS.0.logFC";"R.2-HS.0.logFC";"R.3-HS.0.logFC";"R.4-HS.0.logFC";
    //        "HS.2-HS.1.logFC";"HS.3-HS.2.logFC";"R.1-HS.3.logFC";"R.2-R.1.logFC" 
                                                                                |])>]      [<ChlamyProteome.DoubleArrayConverter>]        Features        : float []
        }

    let csvPathT = sprintf @"%sdata\TranscriptomeData-TuKl\resultTable_hsre.tsv" General.pathToData

    let readerT    = new SchemaReader.Csv.CsvReader<ProteinItemReadT>(schemaMode=SchemaReader.Csv.SchemaMode.Fill)

    /// Variable, contains all raw data from protein DataBase file in csvPath
    let dataProteinT = 
        readerT.ReadFile (csvPathT, General.separatorTab, General.hasHeader) 
        |> Seq.distinctBy (fun x -> x.ProteinGroup)
        |> Seq.toList

    /// Data to find mapman bin for proteins
    let dataMM = 
        readerMM.ReadFile (csvPathMM, General.separatorTab, General.hasHeader) 
        |> Seq.filter (fun x -> x.ProteinIdentifier.Length > 10 && x.ProteinIdentifier.Contains ".")
        |> Seq.toList

    let dataWithBin : ChlamyProteome.ProteinItemRead list =
        dataProteinT
        |> List.choose (fun i -> 
            try Some 
                        {ProteinGroup=i.ProteinGroup; 
                        Bin = 
                            (dataMM 
                            |> List.find (fun ii -> (ii.ProteinIdentifier |> String.subString 2 12)=(i.ProteinGroup.[0] |> String.subString 1 12))).BinString; 
                        Features=i.Features} with
                |ex -> None)


module ArabiTranscriptome =

    // read MapMan terms

    type BinConverter() = 
        inherit SchemaReader.Attribute.ConverterAttribute()
        override this.convertToObj = 
            SchemaReader.Converter.Single (fun (strs : string) -> 
                                                (strs |> String.split ';' |> Array.map (String.replace "GMM:" "" >> String.split '.' ) |> Array.item 0) |> box)

    type IdentifierConverter() = 
        inherit SchemaReader.Attribute.ConverterAttribute()
        override this.convertToObj = 
            SchemaReader.Converter.Single (fun (strs : string) -> 
                                                (strs |> String.map (Char.ToUpper)) |> box)

    type MapManRead = {
        [<SchemaReader.Attribute.FieldAttribute("BINCODE")>]    [<BinConverter>]       BinString           : string []
        [<SchemaReader.Attribute.FieldAttribute("IDENTIFIER")>] [<IdentifierConverter>]       ProteinIdentifier   : string
        }

    let pathMM = sprintf @"%sdata\AT-MapManList.txt" General.pathToData
    
    let readerMM  = new SchemaReader.Csv.CsvReader<MapManRead>(schemaMode=SchemaReader.Csv.SchemaMode.Fill)

    /// Data to find mapman bin for proteins
    let dataMM = 
        readerMM.ReadFile (pathMM, General.separatorTab, General.hasHeader) 
        |> Seq.toArray

    // read csv-file

    type ArabidopsisName() = 
        inherit SchemaReader.Attribute.ConverterAttribute()
        override this.convertToObj = 
            SchemaReader.Converter.Single (fun (strs : string) -> strs |> String.filter (fun c -> c<>'\"') |> box)

    type TripleArrayConverter() = 
        inherit SchemaReader.Attribute.ConverterAttribute()
        override this.convertToObj = 
            SchemaReader.Converter.Collection (fun (strs : seq<string>) -> 
                                                    strs 
                                                    |> Seq.map (fun s -> FSharpAux.String.tryParseFloatDefault nan s) 
                                                    |> Seq.toArray
                                                    |> Array.chunkBySize 3
                                                    //|> Array.map (fun replics -> replics |> Array.filter (fun x -> x<>nan) |> Array.average)
                                                    |> box)

    //// Generate FieldAttribute for AProteinItemReadT
    //let colList =
    //    let exp = "Hlig"
    //    let times = [1;15;180;2880;5760]
    //    let rep = [1;2;3]
    //    let onoff = ["TRUE";"FALSE"]
    //    onoff 
    //    |> List.map (fun d -> 
    //        times 
    //        |> List.map (fun t -> 
    //            rep 
    //            |> List.map (fun r -> sprintf "\\\"%s_t_%i_r_%i_d_%s\\\"" exp t r d))
    //        |> List.concat)
    //    |> List.concat
    //    |> List.toArray
    
    type AProteinItemReadT = {
        [<SchemaReader.Attribute.FieldAttribute("\"\"")>]       [<ArabidopsisName>]           ProteinGroup    : string
        [<SchemaReader.Attribute.FieldAttribute(
        
            //[|"\"Hlig_t_1_r_1_d_TRUE\""; "\"Hlig_t_1_r_2_d_TRUE\"";
            //"\"Hlig_t_1_r_3_d_TRUE\""; "\"Hlig_t_15_r_1_d_TRUE\"";
            //"\"Hlig_t_15_r_2_d_TRUE\""; "\"Hlig_t_15_r_3_d_TRUE\"";
            //"\"Hlig_t_180_r_1_d_TRUE\""; "\"Hlig_t_180_r_2_d_TRUE\"";
            //"\"Hlig_t_180_r_3_d_TRUE\""; "\"Hlig_t_2880_r_1_d_TRUE\"";
            //"\"Hlig_t_2880_r_2_d_TRUE\""; "\"Hlig_t_2880_r_3_d_TRUE\"";
            //"\"Hlig_t_5760_r_1_d_TRUE\""; "\"Hlig_t_5760_r_2_d_TRUE\"";
            //"\"Hlig_t_5760_r_3_d_TRUE\""; "\"Hlig_t_1_r_1_d_FALSE\"";
            //"\"Hlig_t_1_r_2_d_FALSE\""; "\"Hlig_t_1_r_3_d_FALSE\"";
            //"\"Hlig_t_15_r_1_d_FALSE\""; "\"Hlig_t_15_r_2_d_FALSE\"";
            //"\"Hlig_t_15_r_3_d_FALSE\""; "\"Hlig_t_180_r_1_d_FALSE\"";
            //"\"Hlig_t_180_r_2_d_FALSE\""; "\"Hlig_t_180_r_3_d_FALSE\"";
            //"\"Hlig_t_2880_r_1_d_FALSE\""; "\"Hlig_t_2880_r_2_d_FALSE\"";
            //"\"Hlig_t_2880_r_3_d_FALSE\""; "\"Hlig_t_5760_r_1_d_FALSE\"";
            //"\"Hlig_t_5760_r_2_d_FALSE\""; "\"Hlig_t_5760_r_3_d_FALSE\""|])>]      [<TripleArrayConverter>]        Features        : float [] []

            //[|"\"Cold_t_1_r_1_d_TRUE\""; "\"Cold_t_1_r_2_d_TRUE\"";
            //"\"Cold_t_1_r_3_d_TRUE\""; "\"Cold_t_15_r_1_d_TRUE\"";
            //"\"Cold_t_15_r_2_d_TRUE\""; "\"Cold_t_15_r_3_d_TRUE\"";
            //"\"Cold_t_180_r_1_d_TRUE\""; "\"Cold_t_180_r_2_d_TRUE\"";
            //"\"Cold_t_180_r_3_d_TRUE\""; "\"Cold_t_2880_r_1_d_TRUE\"";
            //"\"Cold_t_2880_r_2_d_TRUE\""; "\"Cold_t_2880_r_3_d_TRUE\"";
            //"\"Cold_t_5760_r_1_d_TRUE\""; "\"Cold_t_5760_r_2_d_TRUE\"";
            //"\"Cold_t_5760_r_3_d_TRUE\""; "\"Cold_t_1_r_1_d_FALSE\"";
            //"\"Cold_t_1_r_2_d_FALSE\""; "\"Cold_t_1_r_3_d_FALSE\"";
            //"\"Cold_t_15_r_1_d_FALSE\""; "\"Cold_t_15_r_2_d_FALSE\"";
            //"\"Cold_t_15_r_3_d_FALSE\""; "\"Cold_t_180_r_1_d_FALSE\"";
            //"\"Cold_t_180_r_2_d_FALSE\""; "\"Cold_t_180_r_3_d_FALSE\"";
            //"\"Cold_t_2880_r_1_d_FALSE\""; "\"Cold_t_2880_r_2_d_FALSE\"";
            //"\"Cold_t_2880_r_3_d_FALSE\""; "\"Cold_t_5760_r_1_d_FALSE\"";
            //"\"Cold_t_5760_r_2_d_FALSE\""; "\"Cold_t_5760_r_3_d_FALSE\""|])>]      [<TripleArrayConverter>]        Features        : float [] []
        
            [|"\"Heat_t_1_r_1_d_TRUE\""; "\"Heat_t_1_r_2_d_TRUE\"";
            "\"Heat_t_1_r_3_d_TRUE\""; "\"Heat_t_15_r_1_d_TRUE\"";
            "\"Heat_t_15_r_2_d_TRUE\""; "\"Heat_t_15_r_3_d_TRUE\"";
            "\"Heat_t_180_r_1_d_TRUE\""; "\"Heat_t_180_r_2_d_TRUE\"";
            "\"Heat_t_180_r_3_d_TRUE\""; "\"Heat_t_2880_r_1_d_TRUE\"";
            "\"Heat_t_2880_r_2_d_TRUE\""; "\"Heat_t_2880_r_3_d_TRUE\"";
            "\"Heat_t_5760_r_1_d_TRUE\""; "\"Heat_t_5760_r_2_d_TRUE\"";
            "\"Heat_t_5760_r_3_d_TRUE\""; "\"Heat_t_1_r_1_d_FALSE\"";
            "\"Heat_t_1_r_2_d_FALSE\""; "\"Heat_t_1_r_3_d_FALSE\"";
            "\"Heat_t_15_r_1_d_FALSE\""; "\"Heat_t_15_r_2_d_FALSE\"";
            "\"Heat_t_15_r_3_d_FALSE\""; "\"Heat_t_180_r_1_d_FALSE\"";
            "\"Heat_t_180_r_2_d_FALSE\""; "\"Heat_t_180_r_3_d_FALSE\"";
            "\"Heat_t_2880_r_1_d_FALSE\""; "\"Heat_t_2880_r_2_d_FALSE\"";
            "\"Heat_t_2880_r_3_d_FALSE\""; "\"Heat_t_5760_r_1_d_FALSE\"";
            "\"Heat_t_5760_r_2_d_FALSE\""; "\"Heat_t_5760_r_3_d_FALSE\""|])>]      [<TripleArrayConverter>]        Features        : float [] []
        }

    let csvPathTA = sprintf @"%sdata\Kinetiks\CountData.csv" General.pathToData

    let readerTA    = new SchemaReader.Csv.CsvReader<AProteinItemReadT>(schemaMode=SchemaReader.Csv.SchemaMode.Fill)

    /// Variable, contains all raw data from protein DataBase file in csvPath
    let dataProteinTA = 
        readerTA.ReadFile (csvPathTA, General.separatorComma, General.hasHeader) 
        |> Seq.distinctBy (fun x -> x.ProteinGroup)
        |> Seq.toArray

    // map to MapMan

    //open BioFSharp.BioDB.FaToolDb  

    let mapItemsToMapMan id (item: AProteinItemReadT) = 
        //printfn "name %s" item.ProteinGroup
        let binMM = 
                    /// Calling Web.Server SQL DB
            //(Queries.getOntologyterms MapManOntology item.ProteinGroup).OntologyGroups 
            //|> (fun list -> 
            //                if list.IsEmpty then 
            //                    []
            //                else 
            //                    list
            //                    |> List.head |> fst |> String.split ':' |> Array.item 1 |> String.split '.' |> Array.toList)
            try 
                (dataMM
                |> Array.find (fun i -> i.ProteinIdentifier=item.ProteinGroup)).BinString
            with _ -> [||]
        {
        ID = id;
        ProteinL = [|item.ProteinGroup|];
        OriginalBin = binMM;
        BinL = binMM;
        dataL = 
            item.Features 
            |> Array.map (fun triple -> triple |> Array.choose (fun x -> if x.Equals Double.NaN then None else Some x) |> Array.average)
        } 

    let ItemsWithMapMan = 
        dataProteinTA |> Array.mapi (fun id x -> mapItemsToMapMan id x)
    let itemsWithMapManFound = 
        ItemsWithMapMan |> Array.filter (fun x -> x.OriginalBin<>[||]) |> Array.distinctBy (fun x -> x.ProteinL)
    let itemsWithMapManIdentified =
        itemsWithMapManFound |> Array.filter (fun x -> x.OriginalBin.[0]<>"35")

    /// Filter significant genes

    let pathSign = sprintf @"%sdata\Kinetiks\DESeq2\" General.pathToData 

    let expModeHeat = "Hitze", "Heat"
    let expModeCold = "Kalt", "Cold"
    let expModeLight = "HL", "Hlig"

    let filesSign mode =
        seq
            ["Deseq2-time_8640Controlvs_";
            "Deseq2-time_5940Controlvs_";
            "Deseq2-time_5760Controlvs_";
            "Deseq2-time_2880Controlvs_";
            "Deseq2-time_180Controlvs_";
            "Deseq2-time_11520Controlvs_";
            "Deseq2-time5760vstime_5775_";
            "Deseq2-time5760vstime_5761_";
            "Deseq2-time0vstime_1_";
            "Deseq2-time0vstime_15_"]
        |> Seq.map (fun i -> sprintf "%s%s\%s%s.csv" pathSign (fst mode) i (snd mode))
        
    type DoubleConverter() =
        inherit SchemaReader.Attribute.ConverterAttribute()
        override this.convertToObj =
            SchemaReader.Converter.Single (fun (str : string) -> str |> FSharpAux.String.tryParseFloatDefault nan |> box)

    type AProteinSignT = {
        [<SchemaReader.Attribute.FieldAttribute("\"\"")>]       [<ArabidopsisName>]     ProteinGroup    : string
        [<SchemaReader.Attribute.FieldAttribute("\"pvalue\"")>] [<DoubleConverter>]     Pvalue        : float
        }

    let readerTAsign    = new SchemaReader.Csv.CsvReader<AProteinSignT>(schemaMode=SchemaReader.Csv.SchemaMode.Fill)

    /// Variable, contains all raw data from protein DataBase file in csvPath
    let dataSignFn mode p =
        filesSign mode
        |> Seq.map (fun path ->
            readerTAsign.ReadFile (path, General.separatorComma, General.hasHeader)
            |>  Seq.filter (fun x -> x.Pvalue<=p))
        |> Seq.concat
        |> Seq.distinctBy (fun x -> x.ProteinGroup)
        |> Seq.toArray

    let dataSign = dataSignFn expModeHeat 1. // change here for experiment type and sign level (standard = 0.05, no filter = 1.)

    let filteredAra =
        itemsWithMapManIdentified
        |> Array.filter (fun x -> dataSign |> Array.exists (fun signif -> signif.ProteinGroup=x.ProteinL.[0]))


module ArabiProteome =
    
    let csvPathTA = sprintf @"%sdata\Proteome_201902\SFBcore_Heat_Protein.txt" General.pathToData
    //let csvPathTA = sprintf @"%sdata\Proteome_201902\HeatShock_customNorm.txt" General.pathToData

    type NameConverter() = 
        inherit SchemaReader.Attribute.ConverterAttribute()
        override this.convertToObj = 
            SchemaReader.Converter.Single (fun (strs : string) -> 
                                                (strs |> String.split ';' |> Array.map (fun i -> i |> String.replace "\"" "" |> String.subString 0 9 )) |> box)

    type AProteinItemReadT = {
        [<SchemaReader.Attribute.FieldAttribute("Timepoint")>]       [<NameConverter>]           ProteinGroup    : string []
        [<SchemaReader.Attribute.FieldAttribute(
            [| "A15";"A180";"A2880";"A5760";"DA15";"DA180";"DA2880";"DA5760" |])>]  [<ChlamyProteome.DoubleArrayConverter>] Features    : float [] 
            //[|"0m";"15m";"180m";"2880m";"5760m";"15m";"180m";"2880m";"5760m"|])>]  [<ChlamyProteome.DoubleArrayConverter>] Features    : float [] 
        }

    let readerTA = new SchemaReader.Csv.CsvReader<AProteinItemReadT>(schemaMode=SchemaReader.Csv.SchemaMode.Fill)

    /// Variable, contains all raw data from protein DataBase file in csvPath
    let dataProteinTA = 
        readerTA.ReadFile (csvPathTA, General.separatorTab, General.hasHeader) 
        |> Seq.distinctBy (fun x -> x.ProteinGroup)
        |> Seq.toArray

    let mapItemsToMapMan id (item: AProteinItemReadT) = 
        //printfn "name %s" item.ProteinGroup
        let binMM = 
            [|for a in item.ProteinGroup ->
                try 
                    (ArabiTranscriptome.dataMM |> Array.find (fun i -> i.ProteinIdentifier=a)).BinString 
                with _ -> [||]
            |] 
            //|> Array.filter (fun x -> x.Length>0 && x.[0]<>"35")
            |> fun x ->
                try 
                    (x |> Array.find (fun i -> i<>[||]) )
                with _ -> [||]
        {
        ID = id;
        ProteinL = item.ProteinGroup;
        OriginalBin = binMM;
        BinL = binMM;
        dataL = General.zScoreTransform item.Features 
        } 

    let ItemsWithMapMan = 
        dataProteinTA |> Array.mapi (fun id x -> mapItemsToMapMan id x)

    let itemsWithMapManFound = 
        ItemsWithMapMan |> Array.filter (fun x -> x.OriginalBin<>[||]) |> Array.distinctBy (fun x -> x.ProteinL)

    let itemsWithMapManIdentified =
        itemsWithMapManFound |> Array.filter (fun x -> x.OriginalBin.[0]<>"35")


//let data : Item [] = 
//    (ChlamyProteome.dataAll).[0 .. 40]
//    |> Array.map (fun x -> {x with BinL=[|"0"|]; OriginalBin=[|"0"|]})

//let data : Item [] = 
//    (ArabiTranscriptome.filteredAra).[0 .. 40]
//    |> Array.map (fun x -> {x with BinL=[|"0"|]; OriginalBin=[|"0"|]})
        
    