import pandas as pd
from rdkit import Chem
from rdkit.Chem import MolFromSmiles
from PyQt5.QtWidgets import QMessageBox

class CalcLimit:
    def __init__(self):
        pass

    def validar_nitrosamina(self, smiles):
        try:
            mol = Chem.CanonSmiles(smiles)
            mol = Chem.MolFromSmiles(mol)
            mol = Chem.AddHs(mol)
            pattern = Chem.MolFromSmarts("[#7]-[#7]=[#8]")
            substructure_match = mol.GetSubstructMatch(pattern)
            return substructure_match
        except Exception as e:
            self.show_error_message(e)

    def LerLimite(self, smiles):
        try:
            dictEMA = [{'Name': 'N-Nitroso-Orphenadrine',
            'Smiles': 'CN(CCOC(c1ccccc1)c2ccccc2C)N=O',
            'CPC': '1',
            'AI / ngday-': '18'},
            {'Name': 'N-Nitroso-betahistine',
            'Smiles': 'CN(CCC1=CC=CC=N1)N=O',
            'CPC': '1',
            'AI / ngday-': '18'},
            {'Name': 'N-nitroso-desmethyl-chloropyramine',
            'Smiles': 'CN(CCN(Cc1ccc(Cl)cc1)c2ccccn2)N=O',
            'CPC': '1',
            'AI / ngday-': '18'},
            {'Name': 'N-nitroso-desmethyl-tripelennamine',
            'Smiles': 'CN(CCN(Cc1ccccc1)c2ccccn2)N=O',
            'CPC': '1',
            'AI / ngday-': '18'},
            {'Name': 'P-chlorobenzylamino-pyridine',
            'Smiles': 'Clc1ccc(CN(N=O)c2ccccn2)cc1',
            'CPC': '2',
            'AI / ngday-': '100'},
            {'Name': 'N-nitroso-phenylephrine',
            'Smiles': 'CN(C[C](O)c1cccc(O)c1)N=O',
            'CPC': '2',
            'AI / ngday-': '100'},
            {'Name': 'N-nitroso-rasagiline',
            'Smiles': 'O=NN(CC#C)[C]1CCc2ccccc12',
            'CPC': '2',
            'AI / ngday-': '100'},
            {'Name': 'N-nitroso-sertraline',
            'Smiles': 'CN(N=O)[C]1CC[C](c2ccc(Cl)c(Cl)c2)c3ccccc13',
            'CPC': '2',
            'AI / ngday-': '100'},
            {'Name': '1-cyclopropylmethyl-4-nitrosopiperazine',
            'Smiles': 'O=NN1CCN(CC1)CC2CC2',
            'CPC': '3',
            'AI / ngday-': '400'},
            {'Name': '1-methyl-4-nitroso piperazine',
            'Smiles': 'CN1CCN(CC1)N=O',
            'CPC': '3',
            'AI / ngday-': '400'},
            {'Name': 'nitroso impurity C” [N- (2,6-dimethylphenyl)-2- (4-nitrosopiperazin-1- yl)acetamide]',
            'Smiles': 'Cc1cccc(C)c1NC(=O)CN2CCN(CC2)N=O',
            'CPC': '3',
            'AI / ngday-': '400'},
            {'Name': 'N-nitroso-ambroxol',
            'Smiles': 'C1CC(CCC1N(CC2=C(C(=CC(=C2)Br)Br)N)N=O)O',
            'CPC': '3',
            'AI / ngday-': '400'},
            {'Name': 'N-nitroso-aryl piperazine /N-nitroso-desalkyl quetiapine',
            'Smiles': 'O=NN1CCN(CC1)C2=Nc3ccccc3Sc4ccccc24',
            'CPC': '3',
            'AI / ngday-': '400'},
            {'Name': 'N-nitroso-dabigatran',
            'Smiles': 'Cn1c(CN(N=O)c2ccc(cc2)C(N)=N)nc3cc(ccc13)C(=O)N(CCC(O)=O)c4ccccn4',
            'CPC': '3',
            'AI / ngday-': '400'},
            {'Name': 'N-nitroso-desloratadine',
            'Smiles': 'C1CC2=C(C=CC(=C2)Cl)C(=C3CCN(CC3)N=O)C4=C1C=CC=N4',
            'CPC': '3',
            'AI / ngday-': '400'},
            {'Name': 'N-nitroso-Landiolol',
            'Smiles': 'CC1(C)OC[C](COC(=O)CCc2ccc(OC[C](O)CN(CCNC(=O)N3CCOCC3)N=O)cc2)O1',
            'CPC': '3',
            'AI / ngday-': '400'},
            {'Name': 'N-nitroso-mirabegron',
            'Smiles': 'Nc1scc(CC(=O)Nc2ccc(CCN(C[C](O)c3ccccc3)N=O)cc2)n1',
            'CPC': '3',
            'AI / ngday-': '400'},
            {'Name': 'N-nitroso-N-ethyl-valacyclovir',
            'Smiles': 'CCN(N=O)[C](C(C)C)C(=O)OCCOCn1cnc2C(=O)NC(=Nc12)N',
            'CPC': '3',
            'AI / ngday-': '400'},
            {'Name': 'N-nitroso-N-methyl-valacyclovir',
            'Smiles': 'CC(C)[C](N(C)N=O)C(=O)OCCOCn1cnc2C(=O)NC(=Nc12)N',
            'CPC': '3',
            'AI / ngday-': '400'},
            {'Name': 'N-nitroso-piperazine',
            'Smiles': 'C1CN(CCN1)N=O',
            'CPC': '3',
            'AI / ngday-': '400'},
            {'Name': 'N-nitroso-pramipexole',
            'Smiles': 'CCCN(N=O)[C]1CCc2nc(N)sc2C1',
            'CPC': '3',
            'AI / ngday-': '400'},
            {'Name': 'N-nitroso-trimetazidine',
            'Smiles': 'COC1=C(C(=C(C=C1)CN2CCN(CC2)N=O)OC)OC',
            'CPC': '3',
            'AI / ngday-': '400'},
            {'Name': 'N-nitroso-vortioxetine',
            'Smiles': 'Cc1ccc(Sc2ccccc2N3CCN(CC3)N=O)c(C)c1',
            'CPC': '3',
            'AI / ngday-': '400'},
            {'Name': '1-nitroso-pyrrolopiperidine',
            'Smiles': 'O=NN1CCCC2CNCC12',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'Nitroso-praziquanamine',
            'Smiles': 'O=NN1CC2N(CCc3ccccc23)C(=O)C1',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-N-(2,6-dimethylphenyl)piperidine-2-carboxamide',
            'Smiles': 'Cc1cccc(C)c1NC(=O)C2CCCCN2N=O',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-atenolol',
            'Smiles': 'CC(C)N(CC(COC1=CC=C(C=C1)CC(=O)N)O)N=O',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-bisoprolol',
            'Smiles': 'CC(C)OCCOCc1ccc(OCC(O)CN(N=O)C(C)C)cc1',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-bumetanide',
            'Smiles': 'CCCCN(C1=C(C(=CC(=C1)C(=O)O)S(=O)(=O)N)OC2=CC=CC=C2)N=O',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-ciprofloxacin',
            'Smiles': 'C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCN(CC4)N=O)F)C(=O)O',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-folic acid',
            'Smiles': 'C1=CC(=CC=C1C(=O)NC(CCC(=O)O)C(=O)O)N(CC2=CN=C3C(=N2)C(=O)NC(=N3)N)N=O',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-labetalol',
            'Smiles': 'CC(CCc1ccccc1)N(CC(O)c2ccc(O)c(c2)C(N)=O)N=O',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-levofloxacin',
            'Smiles': 'CC1COC2=C3N1C=C(C(=O)C3=CC(=C2N4CCN(CC4)N=O)F)C(=O)O',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-metoprolol',
            'Smiles': 'CC(C)N(CC(COC1=CC=C(C=C1)CCOC)O)N=O',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-moxifloxacin',
            'Smiles': 'COc1c(N2C[C@@H]3CCCN(N=O)[C@@H]3C2)c(F)cc4C(=O)C(=CN(C5CC5)c14)C(O)=O',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-nebivolol',
            'Smiles': 'C1CC2=C(C=CC(=C2)F)OC1C(CN(CC(C3CCC4=C(O3)C=CC(=C4)F)O)N=O)O',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso Propranolol',
            'Smiles': 'CC(C)N(CC(COC1=CC=CC2=CC=CC=C21)O)N=O',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-pseudoephedrine',
            'Smiles': 'CC(C(C1=CC=CC=C1)O)N(C)N=O',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-sotalol',
            'Smiles': 'CC(C)N(CC(C1=CC=C(C=C1)NS(=O)(=O)C)O)N=O',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-tamsulosin',
            'Smiles': 'CCOc1ccccc1OCCN(N=O)[C](C)Cc2ccc(OC)c(c2)[S](N)(=O)=O',
            'CPC': '4',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-diclofenac',
            'Smiles': 'C1=CC=C(C(=C1)CC(=O)O)N(C2=C(C=CC=C2Cl)Cl)N=O',
            'CPC': '5',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-benazepril',
            'Smiles': 'CCOC(=O)[C](CCc1ccccc1)N(N=O)[C]2CCc3ccccc3N(CC(O)=O)C2=O',
            'CPC': '5',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-bupropion',
            'Smiles': 'CC(C(=O)C1=CC(=CC=C1)Cl)N(C(C)(C)C)N=O',
            'CPC': '5',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-cilazapril',
            'Smiles': 'CCOC(=O)C(CCC1=CC=CC=C1)N(C2CCCN3CCCC(N3C2=O)C(=O)O)N=O',
            'CPC': '5',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-enalapril',
            'Smiles': 'CCOC(=O)C(CCC1=CC=CC=C1)N(C(C)C(=O)N2CCCC2C(=O)O)N=O',
            'CPC': '5',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-lisinopril',
            'Smiles': 'NCCCC[C](N(N=O)[C](CCc1ccccc1)C(O)=O)C(=O)N2CCC[C]2C(O)=O',
            'CPC': '5',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-perindopril',
            'Smiles': 'CCC[C](N(N=O)[C](C)C(=O)N1[C](C[C@@H]2CCCC[C@H]12)C(O)=O)C(=O)OCC',
            'CPC': '5',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-ramipril',
            'Smiles': 'CCOC(=O)[C](CCc1ccccc1)N(N=O)[C](C)C(=O)N2[C](C[C@@H]3CCC[C@H]23)C(O)=O',
            'CPC': '5',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-desmethyl trimebutine',
            'Smiles': 'CCC(COC(=O)c1cc(OC)c(OC)c(OC)c1)(N(C)N=O)c2ccccc2',
            'CPC': '5',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-ketamine',
            'Smiles': 'CN(C1(CCCCC1=O)C2=CC=CC=C2Cl)N=O',
            'CPC': '5',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-salbutamol',
            'Smiles': 'CC(C)(C)N(CC(C1=CC(=C(C=C1)O)CO)O)N=O',
            'CPC': '5',
            'AI / ngday-': '1500'},
            {'Name': 'N-nitroso-vildagliptin',
            'Smiles': 'OC12CC3CC(C1)CC(C3)(C2)N(CC(=O)N4CCC[C]4C#N)N=O',
            'CPC': '5',
            'AI / ngday-': '1500'}]

            df = pd.DataFrame(dictEMA)         
            result = df[df['Smiles'] == smiles]
            #Se na base de dados
            if len(result) > 0:
                cpc = str(result['CPC'].iloc[0])
                ai = str(result['AI / ngday-'].iloc[0])
                return cpc, ai
            else:     
                return None, None
        except Exception as e:
            self.show_error_message(e)
            
    def calcular_H(self, smiles):
        try:
            mol = Chem.CanonSmiles(smiles)
            mol = Chem.MolFromSmiles(mol)
            mol = Chem.AddHs(mol)
            pattern = Chem.MolFromSmarts("[#7]-[#7]=[#8]")
            # Verificar se o padrão está presente na molécula
            substructure_match = mol.GetSubstructMatch(pattern)
            result_dict = {}

            if substructure_match:
                nitrosamine_atoms = set(substructure_match)
                nitrosamine_atoms = list(nitrosamine_atoms)
                # Verificar os vizinhos dos átomos de nitrogênio da nitrosamina
                for atom_idx in nitrosamine_atoms:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    neighbors = atom.GetNeighbors()

                    for neighbor in neighbors:
                        if neighbor.GetAtomicNum() == 6:  # Verificar se é um átomo de carbono
                            carbon_idx = neighbor.GetIdx()
                            hydrogen_count = len([atom for atom in neighbor.GetNeighbors() if atom.GetAtomicNum() == 1])
                            result_dict[carbon_idx] = hydrogen_count
                lista_resultado = list(result_dict.items())
                lado1 = lista_resultado[0]
                lado2 = lista_resultado[1]
                chave_C1, valorH1 = lado1
                chave_C2, valorH2 = lado2
                return valorH1, valorH2, nitrosamine_atoms
        except Exception as e:
            self.show_error_message(e)

    def calcular_terciario(self, smiles):
        try:
            mol = Chem.CanonSmiles(smiles)
            mol = Chem.MolFromSmiles(mol)
            mol = Chem.AddHs(mol)
            
            terciario = Chem.MolFromSmarts("[CX4&$(C(C)(C)(C))][N](N=O)")
            substructure_match = mol.GetSubstructMatch(terciario)
            return substructure_match
        except Exception as e:
            self.show_error_message(e)

    def alfaH_score(self, smiles):
        try:
            scoreH = 0
            smartsetil = Chem.MolFromSmarts('[CH3][CH2]')
            mol = Chem.CanonSmiles(smiles)
            mol = Chem.MolFromSmiles(mol)
            mol = Chem.AddHs(mol)
            valorH1, valorH2, nitrosamine_atoms = self.calcular_H(smiles)
            indice_nitrogenio_nitrosamina = nitrosamine_atoms[0]
            atoms_grupo_etil = mol.GetSubstructMatch(smartsetil)
            # Grupo etil
            vizinhos_grupo_etil = []
            for atom_idx in atoms_grupo_etil:
                atom = mol.GetAtomWithIdx(atom_idx)
                neighbors = atom.GetNeighbors()
                vizinhos_grupo_etil.extend(neighbors)
            #átomos vizinhos do grupo etil é igual ao índice do átomo de nitrogênio da nitrosamina
            grupo_etil_ligado_ao_nitrogenio = any(atom.GetIdx() == indice_nitrogenio_nitrosamina for atom in vizinhos_grupo_etil)
            #contagem dos H
            if (valorH1 == 0 or valorH2 == 0) and (valorH1 == 2 or valorH2 == 2):
                #*  Verificar se o grupo etil está presente na molécula
                if grupo_etil_ligado_ao_nitrogenio:
                    scoreH += 2
                else:
                    scoreH += 3
            elif (valorH1 == 0 or valorH2 == 0) and (valorH1 == 3 or valorH2 == 3):
                scoreH += 2
            elif (valorH1 == 1 or valorH2 == 1) and (valorH1 == 2 or valorH2 == 2):
                scoreH += 3
            elif (valorH1 == 1 or valorH2 == 1) and (valorH1 == 3 or valorH2 == 3):
                scoreH += 3
            elif (valorH1 == 2 and valorH2 == 2):
                scoreH += 1
            elif (valorH1 == 2 or valorH2 == 2) and (valorH1 == 3 or valorH2 == 3):
                scoreH += 1
            else:
                pass
            return valorH1, valorH2, scoreH, grupo_etil_ligado_ao_nitrogenio
        except Exception as e:
            self.show_error_message(e)

    def show_error_message(self, exception):
        if isinstance(exception, Exception):  # Verifica se a variável 'exception' é uma instância de Exception
            error_message = f"Erro: {exception}"
        else:
            error_message = "Não foi possível ler o texto."
        QMessageBox.information(None, "Erro de Leitura", error_message, QMessageBox.Ok)

    

