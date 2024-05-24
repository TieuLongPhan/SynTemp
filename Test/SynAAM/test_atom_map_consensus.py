import unittest
from pathlib import Path

root_dir = Path(__file__).parents[2]
from SynTemp.SynAAM.atom_map_consensus import AAMConsensus


class TestAAMConsensus(unittest.TestCase):
    def setUp(self):
        self.working_dir = "./"
        self.rdt_jar_path = f"{root_dir}/Data/RDT_2.4.1.jar"
        self.mappers = ["local_mapper", "rxn_mapper", "graphormer", "rdt"]
        self.single_dict = {
            "id": "test_1",
            "reactions": "C=O.N#[C-].[K+].O=C(Cl)c1ccccc1>>N#CC(O)C(=O)c1ccccc1",
        }
        self.double_dicts = [
            {
                "id": "test_1",
                "reactions": "C=O.N#[C-].[K+].O=C(Cl)c1ccccc1>>N#CC(O)C(=O)c1ccccc1",
            },
            {
                "id": "test_2",
                "reactions": "Fc1ccc(Br)cc1C(F)(F)F.N#[C-]~[K+].O.O>>O=C(O)c1ccc(Br)cc1C(F)(F)F",
            },
        ]
        self.single_result = {
            "id": "test_1",
            "reactions": "C=O.N#[C-].[K+].O=C(Cl)c1ccccc1>>N#CC(O)C(=O)c1ccccc1",
            "local_mapper": "[CH2:3]=[O:4].[N:1]#[C-:2].[K+].[O:6]=[C:5](Cl)[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1>>[N:1]#[C:2][CH:3]([OH:4])[C:5](=[O:6])[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1",
            "rxn_mapper": "[CH2:3]=[O:4].[N:1]#[C-:2].[K+].[O:6]=[C:5](Cl)[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1>>[N:1]#[C:2][CH:3]([OH:4])[C:5](=[O:6])[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1",
            "graphormer": "[C-:4]#[N:3].[K+:5].[O:2]=[CH2:1].[cH:10]1[cH:11][cH:12][cH:13][cH:14][c:9]1[C:7]([Cl:8])=[O:6]>>[cH:11]1[cH:10][c:9]([cH:14][cH:13][cH:12]1)[C:7](=[O:6])[CH:1]([C:4]#[N:3])[OH:2]",
            "rdt": "[O:1]=[C:2]([Cl:13])[c:3]1[cH:4][cH:5][cH:6][cH:7][cH:8]1.[O:9]=[CH2:10].[C-:11]#[N:12].[K+:14]>>[cH:8]1[cH:7][cH:6][cH:5][cH:4][c:3]1[C:2]([CH:10]([C:11]#[N:12])[OH:9])=[O:1]",
        }

        self.double_results = [
            {
                "id": "test_1",
                "reactions": "C=O.N#[C-].[K+].O=C(Cl)c1ccccc1>>N#CC(O)C(=O)c1ccccc1",
                "local_mapper": "[CH2:3]=[O:4].[N:1]#[C-:2].[K+].[O:6]=[C:5](Cl)[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1>>[N:1]#[C:2][CH:3]([OH:4])[C:5](=[O:6])[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1",
                "rxn_mapper": "[CH2:3]=[O:4].[N:1]#[C-:2].[K+].[O:6]=[C:5](Cl)[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1>>[N:1]#[C:2][CH:3]([OH:4])[C:5](=[O:6])[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1",
                "graphormer": "[C-:4]#[N:3].[K+:5].[O:2]=[CH2:1].[cH:10]1[cH:11][cH:12][cH:13][cH:14][c:9]1[C:7]([Cl:8])=[O:6]>>[cH:11]1[cH:10][c:9]([cH:14][cH:13][cH:12]1)[C:7](=[O:6])[CH:1]([C:4]#[N:3])[OH:2]",
                "rdt": "[O:1]=[C:2]([Cl:13])[c:3]1[cH:4][cH:5][cH:6][cH:7][cH:8]1.[O:9]=[CH2:10].[C-:11]#[N:12].[K+:14]>>[cH:8]1[cH:7][cH:6][cH:5][cH:4][c:3]1[C:2]([CH:10]([C:11]#[N:12])[OH:9])=[O:1]",
            },
            {
                "id": "test_2",
                "reactions": "Fc1ccc(Br)cc1C(F)(F)F.N#[C-]~[K+].O.O>>O=C(O)c1ccc(Br)cc1C(F)(F)F",
                "local_mapper": "F[c:4]1[cH:5][cH:6][c:7]([Br:8])[cH:9][c:10]1[C:11]([F:12])([F:13])[F:14].N#[C-:2]~[K+].[OH2:1].[OH2:3]>>[O:1]=[C:2]([OH:3])[c:4]1[cH:5][cH:6][c:7]([Br:8])[cH:9][c:10]1[C:11]([F:12])([F:13])[F:14]",
                "rxn_mapper": "F[c:4]1[cH:5][cH:6][c:7]([Br:8])[cH:9][c:10]1[C:11]([F:12])([F:13])[F:14].N#[C-:2]~[K+].[OH2:1].[OH2:3]>>[O:1]=[C:2]([OH:3])[c:4]1[cH:5][cH:6][c:7]([Br:8])[cH:9][c:10]1[C:11]([F:12])([F:13])[F:14]",
                "graphormer": "[F:10][C:9]([F:11])([F:12])[c:8]1[cH:7][c:5]([Br:6])[cH:4][cH:3][c:2]1[F:1].[K+:15]~[C-:14]#[N:13].[OH2:16].[OH2:17]>>[F:10][C:9]([F:11])([F:12])[c:8]1[cH:7][c:5]([Br:6])[cH:4][cH:3][c:2]1[C:14](=[O:17])[OH:16]",
                "rdt": "Fc1ccc(Br)cc1C(F)(F)F.N#[C-]~[K+].O.O>>O=C(O)c1ccc(Br)cc1C(F)(F)F",
            },
        ]  # cannot process

    def test_single_consensus(self):
        aam = AAMConsensus(self.single_dict, self.mappers)
        results = aam.single_consensus(
            self.single_dict, "reactions", self.rdt_jar_path, self.working_dir
        )
        self.assertEqual(results, self.single_result)

    def test_batch_consensus_unsafe_mode(self):
        aam = AAMConsensus(self.double_dicts, self.mappers)
        results = aam.batch_consensus(
            self.double_dicts,
            "reactions",
            2,
            10,
            False,
            self.rdt_jar_path,
            self.working_dir,
        )
        self.assertEqual(results, self.double_results)

    def test_batch_consensus_safe_mode(self):
        aam = AAMConsensus(self.double_dicts, self.mappers)
        long_smiles = [
            {
                "id": "test_3",
                "reactions": "NC1C2=CC(=C(OCC(O)COC(=O)C3=C(C4=CC=CC=C4)C4C=CC3C4)C=C2OCC(O)COC(=O)C2=C(C3=CC=CC=C3)C3C=CC2C3)C(C)C2=CC(=C(OCC(O)COC(=O)C3=C(C4=CC=CC=C4)C4C=CC3C4)C=C2OCC(O)COC(=O)C2=C(C3=CC=CC=C3)C3C=CC2C3)C(C)C2=CC(=C(OCC(O)COC(=O)C3=C(C4=CC=CC=C4)C4C=CC3C4)C=C2OCC(O)COC(=O)C2=C(C3=CC=CC=C3)C3C=CC2C3)C(C)C2=CC1=C(OCC(O)COC(=O)C1=C(C3=CC=CC=C3)C3C=CC1C3)C=C2OCC(O)COC(=O)C1=C(C2=CC=CC=C2)C2C=CC1C2>>NC1C2=CC(=C(OCC(O)COC(=O)C34C5CC6C(C53)C64C3=CC=CC=C3)C=C2OCC(O)COC(=O)C23C4CC5C(C42)C53C2=CC=CC=C2)C(C)C2=CC(=C(OCC(O)COC(=O)C34C5CC6C(C53)C64C3=CC=CC=C3)C=C2OCC(O)COC(=O)C23C4CC5C(C42)C53C2=CC=CC=C2)C(C)C2=CC(=C(OCC(O)COC(=O)C34C5CC6C(C53)C64C3=CC=CC=C3)C=C2OCC(O)COC(=O)C23C4CC5C(C42)C53C2=CC=CC=C2)C(C)C2=CC1=C(OCC(O)COC(=O)C13C4CC5C(C41)C53C1=CC=CC=C1)C=C2OCC(O)COC(=O)C12C3CC4C(C31)C42C1=CC=CC=C1",
            }
        ]
        results = aam.batch_consensus(
            long_smiles, "reactions", 1, 10, True, self.rdt_jar_path, self.working_dir
        )
        expected_list = [
            {
                "id": "test_3",
                "reactions": "NC1C2=CC(=C(OCC(O)COC(=O)C3=C(C4=CC=CC=C4)C4C=CC3C4)C=C2OCC(O)COC(=O)C2=C(C3=CC=CC=C3)C3C=CC2C3)C(C)C2=CC(=C(OCC(O)COC(=O)C3=C(C4=CC=CC=C4)C4C=CC3C4)C=C2OCC(O)COC(=O)C2=C(C3=CC=CC=C3)C3C=CC2C3)C(C)C2=CC(=C(OCC(O)COC(=O)C3=C(C4=CC=CC=C4)C4C=CC3C4)C=C2OCC(O)COC(=O)C2=C(C3=CC=CC=C3)C3C=CC2C3)C(C)C2=CC1=C(OCC(O)COC(=O)C1=C(C3=CC=CC=C3)C3C=CC1C3)C=C2OCC(O)COC(=O)C1=C(C2=CC=CC=C2)C2C=CC1C2>>NC1C2=CC(=C(OCC(O)COC(=O)C34C5CC6C(C53)C64C3=CC=CC=C3)C=C2OCC(O)COC(=O)C23C4CC5C(C42)C53C2=CC=CC=C2)C(C)C2=CC(=C(OCC(O)COC(=O)C34C5CC6C(C53)C64C3=CC=CC=C3)C=C2OCC(O)COC(=O)C23C4CC5C(C42)C53C2=CC=CC=C2)C(C)C2=CC(=C(OCC(O)COC(=O)C34C5CC6C(C53)C64C3=CC=CC=C3)C=C2OCC(O)COC(=O)C23C4CC5C(C42)C53C2=CC=CC=C2)C(C)C2=CC1=C(OCC(O)COC(=O)C13C4CC5C(C41)C53C1=CC=CC=C1)C=C2OCC(O)COC(=O)C12C3CC4C(C31)C42C1=CC=CC=C1",
                "local_mapper": "NC1C2=CC(=C(OCC(O)COC(=O)C3=C(C4=CC=CC=C4)C4C=CC3C4)C=C2OCC(O)COC(=O)C2=C(C3=CC=CC=C3)C3C=CC2C3)C(C)C2=CC(=C(OCC(O)COC(=O)C3=C(C4=CC=CC=C4)C4C=CC3C4)C=C2OCC(O)COC(=O)C2=C(C3=CC=CC=C3)C3C=CC2C3)C(C)C2=CC(=C(OCC(O)COC(=O)C3=C(C4=CC=CC=C4)C4C=CC3C4)C=C2OCC(O)COC(=O)C2=C(C3=CC=CC=C3)C3C=CC2C3)C(C)C2=CC1=C(OCC(O)COC(=O)C1=C(C3=CC=CC=C3)C3C=CC1C3)C=C2OCC(O)COC(=O)C1=C(C2=CC=CC=C2)C2C=CC1C2>>NC1C2=CC(=C(OCC(O)COC(=O)C34C5CC6C(C53)C64C3=CC=CC=C3)C=C2OCC(O)COC(=O)C23C4CC5C(C42)C53C2=CC=CC=C2)C(C)C2=CC(=C(OCC(O)COC(=O)C34C5CC6C(C53)C64C3=CC=CC=C3)C=C2OCC(O)COC(=O)C23C4CC5C(C42)C53C2=CC=CC=C2)C(C)C2=CC(=C(OCC(O)COC(=O)C34C5CC6C(C53)C64C3=CC=CC=C3)C=C2OCC(O)COC(=O)C23C4CC5C(C42)C53C2=CC=CC=C2)C(C)C2=CC1=C(OCC(O)COC(=O)C13C4CC5C(C41)C53C1=CC=CC=C1)C=C2OCC(O)COC(=O)C12C3CC4C(C31)C42C1=CC=CC=C1",
                "rxn_mapper": "NC1C2=CC(=C(OCC(O)COC(=O)C3=C(C4=CC=CC=C4)C4C=CC3C4)C=C2OCC(O)COC(=O)C2=C(C3=CC=CC=C3)C3C=CC2C3)C(C)C2=CC(=C(OCC(O)COC(=O)C3=C(C4=CC=CC=C4)C4C=CC3C4)C=C2OCC(O)COC(=O)C2=C(C3=CC=CC=C3)C3C=CC2C3)C(C)C2=CC(=C(OCC(O)COC(=O)C3=C(C4=CC=CC=C4)C4C=CC3C4)C=C2OCC(O)COC(=O)C2=C(C3=CC=CC=C3)C3C=CC2C3)C(C)C2=CC1=C(OCC(O)COC(=O)C1=C(C3=CC=CC=C3)C3C=CC1C3)C=C2OCC(O)COC(=O)C1=C(C2=CC=CC=C2)C2C=CC1C2>>NC1C2=CC(=C(OCC(O)COC(=O)C34C5CC6C(C53)C64C3=CC=CC=C3)C=C2OCC(O)COC(=O)C23C4CC5C(C42)C53C2=CC=CC=C2)C(C)C2=CC(=C(OCC(O)COC(=O)C34C5CC6C(C53)C64C3=CC=CC=C3)C=C2OCC(O)COC(=O)C23C4CC5C(C42)C53C2=CC=CC=C2)C(C)C2=CC(=C(OCC(O)COC(=O)C34C5CC6C(C53)C64C3=CC=CC=C3)C=C2OCC(O)COC(=O)C23C4CC5C(C42)C53C2=CC=CC=C2)C(C)C2=CC1=C(OCC(O)COC(=O)C13C4CC5C(C41)C53C1=CC=CC=C1)C=C2OCC(O)COC(=O)C12C3CC4C(C31)C42C1=CC=CC=C1",
                "graphormer": "[O:135]([C:136](=[O:137])[C:138]1=[C:139]([CH:146]2[CH2:150][CH:149]1[CH:148]=[CH:147]2)[C:140]=1[CH:145]=[CH:144][CH:143]=[CH:142][CH:141]=1)[CH2:134][CH:132]([CH2:131][O:130][C:129]1=[CH:128][C:106](=[C:105]2[CH:104]=[C:103]1[CH:101]([C:55]=1[CH:54]=[C:53]([CH:51]([C:5]3=[C:6]([O:7][CH2:8][CH:9]([OH:10])[CH2:11][O:12][C:13](=[O:14])[C:15]4=[C:16]([CH:23]5[CH:24]=[CH:25][CH:26]4[CH2:27]5)[C:17]=4[CH:22]=[CH:21][CH:20]=[CH:19][CH:18]=4)[CH:28]=[C:29]([C:3]([CH:2]([NH2:1])[C:155]4=[C:156]([O:157][CH2:158][CH:159]([CH2:161][O:162][C:163]([C:165]=5[CH:176]6[CH2:177][CH:173]([C:166]=5[C:167]=5[CH:172]=[CH:171][CH:170]=[CH:169][CH:168]=5)[CH:174]=[CH:175]6)=[O:164])[OH:160])[CH:178]=[C:179]([O:180][CH2:181][CH:182]([CH2:184][O:185][C:186]([C:188]=5[CH:199]6[CH:198]=[CH:197][CH:196]([C:189]=5[C:190]5=[CH:191][CH:192]=[CH:193][CH:194]=[CH:195]5)[CH2:200]6)=[O:187])[OH:183])[C:153]([CH:151]2[CH3:152])=[CH:154]4)=[CH:4]3)[O:30][CH2:31][CH:32]([CH2:34][O:35][C:36]([C:38]2=[C:39]([CH:46]3[CH2:50][CH:49]2[CH:48]=[CH:47]3)[C:40]2=[CH:41][CH:42]=[CH:43][CH:44]=[CH:45]2)=[O:37])[OH:33])[CH3:52])[C:79]([O:80][CH2:81][CH:82]([OH:83])[CH2:84][O:85][C:86]([C:88]=2[CH:99]3[CH2:100][CH:96]([CH:97]=[CH:98]3)[C:89]=2[C:90]=2[CH:95]=[CH:94][CH:93]=[CH:92][CH:91]=2)=[O:87])=[CH:78][C:56]=1[O:57][CH2:58][CH:59]([CH2:61][O:62][C:63](=[O:64])[C:65]1=[C:66]([C:67]2=[CH:68][CH:69]=[CH:70][CH:71]=[CH:72]2)[CH:73]2[CH:74]=[CH:75][CH:76]1[CH2:77]2)[OH:60])[CH3:102])[O:107][CH2:108][CH:109]([CH2:111][O:112][C:113]([C:115]=1[CH:126]2[CH:125]=[CH:124][CH:123]([C:116]=1[C:117]=1[CH:122]=[CH:121][CH:120]=[CH:119][CH:118]=1)[CH2:127]2)=[O:114])[OH:110])[OH:133]>>[C:116]12([C:115]3([CH:126]4[CH:127]3[CH:123]1[CH:124]2[CH2:125]4)[C:113](=[O:114])[O:112][CH2:111][CH:109]([CH2:108][O:107][C:106]=1[CH:128]=[C:129]([C:103]=2[CH:101]([C:55]=3[CH:54]=[C:53]([C:79]([O:80][CH2:81][CH:82]([CH2:84][O:85][C:86]([C:88]45[C:89]6([CH:96]7[CH:97]6[CH2:98][CH:99]4[CH:100]57)[C:90]=4[CH:95]=[CH:94][CH:93]=[CH:92][CH:91]=4)=[O:87])[OH:83])=[CH:78][C:56]=3[O:57][CH2:58][CH:59]([CH2:61][O:62][C:63]([C:65]34[CH:76]5[CH2:75][CH:74]6[C:66]3([C:67]=3[CH:68]=[CH:69][CH:70]=[CH:71][CH:72]=3)[CH:73]6[CH:77]45)=[O:64])[OH:60])[CH:51]([CH3:52])[C:5]3=[C:6]([O:7][CH2:8][CH:9]([OH:10])[CH2:11][O:12][C:13]([C:15]45[CH:26]6[CH2:25][CH:24]7[CH:23]([C:16]47[C:17]4=[CH:18][CH:19]=[CH:20][CH:21]=[CH:22]4)[CH:27]56)=[O:14])[CH:28]=[C:29]([C:3](=[CH:4]3)[CH:2]([NH2:1])[C:155]3=[C:156]([O:157][CH2:158][CH:159]([OH:160])[CH2:161][O:162][C:163]([C:165]45[C:166]6([C:167]7=[CH:172][CH:171]=[CH:170][CH:169]=[CH:168]7)[CH:174]7[CH2:175][CH:176]4[CH:177]5[CH:173]67)=[O:164])[CH:178]=[C:179]([C:153](=[CH:154]3)[CH:151]([CH3:152])[C:105]=1[CH:104]=2)[O:180][CH2:181][CH:182]([CH2:184][O:185][C:186]([C:188]12[CH:200]3[CH:196]4[C:189]1([C:190]1=[CH:191][CH:192]=[CH:193][CH:194]=[CH:195]1)[CH:197]4[CH2:198][CH:199]23)=[O:187])[OH:183])[O:30][CH2:31][CH:32]([OH:33])[CH2:34][O:35][C:36]([C:38]12[C:39]3([C:40]4=[CH:41][CH:42]=[CH:43][CH:44]=[CH:45]4)[CH:46]4[CH:50]1[CH:49]2[CH2:48][CH:47]34)=[O:37])[CH3:102])[O:130][CH2:131][CH:132]([OH:133])[CH2:134][O:135][C:136](=[O:137])[C:138]12[C:139]3([CH:146]4[CH:150]1[CH:149]2[CH2:148][CH:147]34)[C:140]1=[CH:141][CH:142]=[CH:143][CH:144]=[CH:145]1)[OH:110])[C:117]1=[CH:118][CH:119]=[CH:120][CH:121]=[CH:122]1",
                "rdt": "[O:1]=[C:2]([O:3][CH2:4][CH:5]([OH:6])[CH2:7][O:8][c:9]1[cH:10][c:11]([O:12][CH2:13][CH:14]([OH:15])[CH2:16][O:17][C:18](=[O:19])[C:20]2=[C:21]([c:22]3[cH:23][cH:24][cH:25][cH:26][cH:27]3)[CH:28]4[CH:29]=[CH:30][CH:31]2[CH2:32]4)[c:33]5[cH:34][c:35]1[CH:36]([NH2:37])[c:38]6[cH:39][c:40]([c:41]([O:42][CH2:43][CH:44]([OH:45])[CH2:46][O:47][C:48](=[O:49])[C:50]7=[C:51]([c:52]8[cH:53][cH:54][cH:55][cH:56][cH:57]8)[CH:58]9[CH:59]=[CH:60][CH:61]7[CH2:62]9)[cH:63][c:64]6[O:65][CH2:66][CH:67]([OH:68])[CH2:69][O:70][C:71](=[O:72])[C:73]%10=[C:74]([c:75]%11[cH:76][cH:77][cH:78][cH:79][cH:80]%11)[CH:81]%12[CH:82]=[CH:83][CH:84]%10[CH2:85]%12)[CH:86]([c:87]%13[cH:88][c:89]([c:90]([O:91][CH2:92][CH:93]([OH:94])[CH2:95][O:96][C:97](=[O:98])[C:99]%14=[C:100]([c:101]%15[cH:102][cH:103][cH:104][cH:105][cH:106]%15)[CH:107]%16[CH:108]=[CH:109][CH:110]%14[CH2:111]%16)[cH:112][c:113]%13[O:114][CH2:115][CH:116]([OH:117])[CH2:118][O:119][C:120](=[O:121])[C:122]%17=[C:123]([c:124]%18[cH:125][cH:126][cH:127][cH:128][cH:129]%18)[CH:130]%19[CH:131]=[CH:132][CH:133]%17[CH2:134]%19)[CH:135]([c:136]%20[cH:137][c:138]([c:139]([O:140][CH2:141][CH:142]([OH:143])[CH2:144][O:145][C:146](=[O:147])[C:148]%21=[C:149]([c:150]%22[cH:151][cH:152][cH:153][cH:154][cH:155]%22)[CH:156]%23[CH:157]=[CH:158][CH:159]%21[CH2:160]%23)[cH:161][c:162]%20[O:163][CH2:164][CH:165]([OH:166])[CH2:167][O:168][C:169](=[O:170])[C:171]%24=[C:172]([c:173]%25[cH:174][cH:175][cH:176][cH:177][cH:178]%25)[CH:179]%26[CH:180]=[CH:181][CH:182]%24[CH2:183]%26)[CH:184]5[CH3:185])[CH3:186])[CH3:187])[C:188]%27=[C:189]([c:190]%28[cH:191][cH:192][cH:193][cH:194][cH:195]%28)[CH:196]%29[CH:197]=[CH:198][CH:199]%27[CH2:200]%29>>[O:1]=[C:2]([O:3][CH2:4][CH:5]([OH:6])[CH2:7][O:8][c:9]1[cH:10][c:11]([O:12][CH2:13][CH:14]([OH:15])[CH2:16][O:17][C:18](=[O:19])[C:20]23[CH:30]4[CH2:29][CH:28]5[CH:32]([CH:31]24)[C:21]35[c:22]6[cH:23][cH:24][cH:25][cH:26][cH:27]6)[c:33]7[cH:34][c:35]1[CH:36]([NH2:37])[c:38]8[cH:39][c:40]([c:41]([O:42][CH2:43][CH:44]([OH:45])[CH2:46][O:47][C:48](=[O:49])[C:50]9%10[CH:60]%11[CH2:59][CH:58]%12[CH:62]([CH:61]9%11)[C:51]%10%12[c:52]%13[cH:53][cH:54][cH:55][cOutput",
            }
        ]
        self.assertEqual(results, expected_list)
