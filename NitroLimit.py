import sys
import os

files_dir = os.path.join(os.path.dirname(__file__), "Files")
sys.path.append(files_dir)

from main_ui import Ui_MainWindow
from functions import CalcLimit
from PyQt5.QtWidgets import QMainWindow, QDialog, QMessageBox, QApplication 
from PyQt5.QtGui import QImage, QPixmap
from io import BytesIO
from rdkit.Chem.Draw import MolToImage
from rdkit import Chem

class App1(QMainWindow, QDialog):
    def __init__(self):
        super(App1, self).__init__()
        self.gui = Ui_MainWindow()
        self.calc = CalcLimit()
        self.gui.setupUi(self)
        #imagens
        acido = self.show_image('C(=O)O')
        self.gui.label_acido.setPixmap(acido)
        pirrolidina = self.show_image('C1CCCN1N=O')
        self.gui.label_pirrol.setPixmap(pirrolidina)
        anel_6S = self.show_image('C1CSCCN1N=O')
        self.gui.label_S.setPixmap(anel_6S)
        nitroso_5ou6 = self.show_image('C1CNCCN1N=O')
        self.gui.label_nitroso56.setPixmap(nitroso_5ou6)
        nitroso_morf = self.show_image('C1COCCN1N=O')
        self.gui.label_morfolina.setPixmap(nitroso_morf)
        nitro7 = self.show_image('C1CCCCCN1N=O')
        self.gui.label_nitro7.setPixmap(nitro7)
        arila = self.show_image('C(C1=CC=CO1)N(N=O)C')
        self.gui.label_arila.setPixmap(arila)
        metila = self.show_image('C1(CCCN(C1)N=O)C')
        self.gui.label_metila.setPixmap(metila)
        more5C = self.show_image('C(CN(N=O)C1=CC=C(O1)CC)C2=CC=CC=C2')
        self.gui.label_5Cmore.setPixmap(more5C)
        retirador = self.show_image('C(CN(C(CC(C)C)[S](C)(=O)=O)N=O)(N)=O')
        self.gui.label_retirador.setPixmap(retirador)
        hidroxila = self.show_image('C(CN(CC(C(C)C)O[H])N=O)(C(C)C)O[H]')
        self.gui.label_hidroxila.setPixmap(hidroxila) 
        
        self.gui.btn_exibir.clicked.connect(self.fCalcularLimites)
        self.gui.btn_calcular.clicked.connect(self.fCalcularPotencia)
        self.gui.btn_clear.clicked.connect(self.fLimparAll)

                   
    def fLerSmiles(self):   ## smiles
        try:
            nm_smiles = self.gui.ln_smiles.text()
            return nm_smiles
        except Exception as e:
            self.show_error_message(e)          
    
    def fLimparCampos(self):
        try:
            self.gui.h_alfa.setText(' ')
            self.gui.mais_halfa.setText(' ')
            self.gui.carbono_alfa.setText(' ')
            self.gui.h_alfa1.setText(' ')
            self.gui.h_alfa2.setText(' ')
            self.gui.score_halfalados.setText(' ')
            self.gui.categoria.setText(' ')
            self.gui.limite_diario.setText(' ')
            self.gui.score_final.setText(' ')
            self.gui.ln_limite_util.setText(' ')
            self.gui.chk_acido.setChecked(False)
            self.gui.chk_pirrolidina.setChecked(False)
            self.gui.chk_anelS.setChecked(False)
            self.gui.chk_arila.setChecked(False)
            self.gui.chk_5ou6.setChecked(False)
            self.gui.chk_morfolina.setChecked(False)
            self.gui.chk_7membros.setChecked(False)
            self.gui.chk_metil.setChecked(False)
            self.gui.chk_5C.setChecked(False)
            self.gui.chk_retirador1.setChecked(False)
            self.gui.chk_retirador2.setChecked(False)
            self.gui.chk_oh1.setChecked(False)
            self.gui.chk_oh2.setChecked(False)
            self.gui.ln_acidocarb.setText(' ')
            self.gui.ln_pirrolidina.setText(' ')
            self.gui.ln_anelS.setText(' ')
            self.gui.ln_arila.setText(' ')
            self.gui.ln_5ou6.setText(' ')
            self.gui.ln_morfolna.setText(' ')
            self.gui.ln_7membros.setText(' ')
            self.gui.ln_metil.setText(' ')
            self.gui.ln_5C.setText(' ')
            self.gui.score_ret.setText(' ')
            self.gui.ln_hidroxila.setText(' ')            
            self.gui.ln_grupoetil.setText(' ')
        except Exception as e:
            self.show_error_message(e)    

    def fLimparAll(self):
        try:
            self.gui.ln_smiles.clear()
            self.gui.label_smiles.setText(' ')
            self.gui.h_alfa.setText(' ')
            self.gui.mais_halfa.setText(' ')
            self.gui.carbono_alfa.setText(' ')
            self.gui.h_alfa1.setText(' ')
            self.gui.h_alfa2.setText(' ')
            self.gui.score_halfalados.setText(' ')
            self.gui.categoria.setText(' ')
            self.gui.limite_diario.setText(' ')
            self.gui.score_final.setText(' ')
            self.gui.ln_limite_util.setText(' ')
            self.gui.chk_acido.setChecked(False)
            self.gui.chk_pirrolidina.setChecked(False)
            self.gui.chk_anelS.setChecked(False)
            self.gui.chk_arila.setChecked(False)
            self.gui.chk_5ou6.setChecked(False)
            self.gui.chk_morfolina.setChecked(False)
            self.gui.chk_7membros.setChecked(False)
            self.gui.chk_metil.setChecked(False)
            self.gui.chk_5C.setChecked(False)
            self.gui.chk_retirador1.setChecked(False)
            self.gui.chk_retirador2.setChecked(False)
            self.gui.chk_oh1.setChecked(False)
            self.gui.chk_oh2.setChecked(False)
            self.gui.ln_acidocarb.setText(' ')
            self.gui.ln_pirrolidina.setText(' ')
            self.gui.ln_anelS.setText(' ')
            self.gui.ln_arila.setText(' ')
            self.gui.ln_5ou6.setText(' ')
            self.gui.ln_morfolna.setText(' ')
            self.gui.ln_7membros.setText(' ')
            self.gui.ln_metil.setText(' ')
            self.gui.ln_5C.setText(' ')
            self.gui.score_ret.setText(' ')
            self.gui.ln_hidroxila.setText(' ')            
            self.gui.ln_grupoetil.setText(' ')
        except Exception as e:
            self.show_error_message(e)    

    def show_image(self, smiles):
        try:
            image = Chem.MolFromSmiles(smiles)
            img = MolToImage(image)
            width, height = 90, 90
            img = img.resize((width, height))
            image_buffer = BytesIO()
            img.save(image_buffer, format='PNG')
            image_buffer.seek(0)
            qimg = QImage.fromData(image_buffer.getvalue())
            pixmap = QPixmap.fromImage(qimg)
            #scaled_pixmap = pixmap.scaled(self.gui.label_2.size(), aspectRatioMode=Qt.KeepAspectRatio, transformMode=Qt.SmoothTransformation)
            return pixmap
        except Exception as e:
            self.show_error_message(e)  

    def show_image_smiles(self, smiles):
        try:
            image = Chem.MolFromSmiles(smiles)
            img = MolToImage(image)
            width, height = 250, 200
            img = img.resize((width, height))
            image_buffer = BytesIO()
            img.save(image_buffer, format='PNG')
            image_buffer.seek(0)
            qimg = QImage.fromData(image_buffer.getvalue())
            pixmap = QPixmap.fromImage(qimg)
            #scaled_pixmap = pixmap.scaled(self.gui.label_2.size(), aspectRatioMode=Qt.KeepAspectRatio, transformMode=Qt.SmoothTransformation)
            return pixmap
        except Exception as e:
            self.show_error_message(e)  

    def fCalcularLimites(self):
        try:
            self.fLimparCampos()
            smiles = self.fLerSmiles()
            substructure_match = self.calc.validar_nitrosamina(smiles)
            if substructure_match:
                scaled_pixmap = self.show_image_smiles(smiles)
                self.gui.label_smiles.setPixmap(scaled_pixmap)
                cpc, ai  = self.calc.LerLimite(smiles)
                if cpc:           
                    self.gui.categoria.setText(str(cpc))
                    self.gui.limite_diario.setText(str(ai))
                    self.gui.ln_limite_util.setText('Limite Conhecido EMA')
                else:
                    self.fCalcularScore5()                   
        except Exception as e:
            self.show_error_message(e)  
        
    def fCalcularScore5(self):
        try:
            smiles = self.fLerSmiles()
            substructure_match = self.calc.validar_nitrosamina(smiles)
            if substructure_match:
                valorH1, valorH2, nitrosamine_atoms = self.calc.calcular_H(smiles)
                if valorH1 > 0 or valorH2 > 0:
                    self.gui.h_alfa.setText("Sim")
                    if (valorH1 > 1 and valorH2 >= 0) or (valorH1 >= 0 and valorH2 > 1):
                        self.gui.mais_halfa.setText("Sim")
                        terciario = self.calc.calcular_terciario(smiles)
                        if terciario:
                            self.gui.carbono_alfa.setText("Sim")
                            self.gui.categoria.setText(str(5))
                            self.gui.limite_diario.setText(str(1500))
                            self.gui.ln_limite_util.setText('Limite Calculado')
                        else:
                            self.gui.carbono_alfa.setText("Não")
                            self.gui.ln_limite_util.setText('Calcule a potência')
                    else:
                        self.gui.mais_halfa.setText("Não")
                        self.gui.categoria.setText(str(5))
                        self.gui.limite_diario.setText(str(1500))
                        self.gui.ln_limite_util.setText('Limite Calculado')
                else:
                    self.gui.h_alfa.setText("Não")
                    self.gui.categoria.setText(str(5))
                    self.gui.limite_diario.setText(str(1500))
                    self.gui.ln_limite_util.setText('Limite Calculado')
        except Exception as e:
            self.show_error_message(e)  

    def fLerCheckbox(self):
        try:
            score = 0
            grupos = []
            #ativante
            chk_acido = lambda x: True if x == 2 else False
            chk_acido = chk_acido(self.gui.chk_acido.checkState())
            if chk_acido:
                grupos.append('chk_acido')
                score += 3
            chk_pirrolidina = lambda x: True if x == 2 else False
            chk_pirrolidina = chk_pirrolidina(self.gui.chk_pirrolidina.checkState())
            if chk_pirrolidina:
                grupos.append('chk_pirrolidina')
                score += 3
            chk_anelS = lambda x: True if x == 2 else False
            chk_anelS = chk_anelS(self.gui.chk_anelS.checkState())
            if chk_anelS:
                grupos.append('chk_anelS')
                score += 3
            chk_5ou6 = lambda x: True if x == 2 else False
            chk_5ou6 = chk_5ou6(self.gui.chk_5ou6.checkState())
            if chk_5ou6:
                grupos.append('chk_5ou6')
                score += 2
            chk_morfolina = lambda x: True if x == 2 else False
            chk_morfolina = chk_morfolina(self.gui.chk_morfolina.checkState())
            if chk_morfolina:
                grupos.append('chk_morfolina')
                score += 1
            chk_7membros = lambda x: True if x == 2 else False
            chk_7membros = chk_7membros(self.gui.chk_7membros.checkState())
            if chk_7membros:
                grupos.append('chk_7membros')
                score += 1
            chk_5C = lambda x: True if x == 2 else False
            chk_5C = chk_5C(self.gui.chk_5C.checkState())
            if chk_5C:
                grupos.append('chk_5C')
                score += 1
            chk_retirador1 = lambda x: True if x == 2 else False
            chk_retirador1 = chk_retirador1(self.gui.chk_retirador1.checkState())
            if chk_retirador1:
                grupos.append('chk_retirador1')
                score += 1
            chk_retirador2 = lambda x: True if x == 2 else False
            chk_retirador2 = chk_retirador2(self.gui.chk_retirador2.checkState())
            if chk_retirador2:
                grupos.append('chk_retirador2')
                score += 2
            chk_oh1 = lambda x: True if x == 2 else False
            chk_oh1 = chk_oh1(self.gui.chk_oh1.checkState())
            if chk_oh1:
                grupos.append('chk_oh1')
                score += 1
            chk_oh2 = lambda x: True if x == 2 else False
            chk_oh2 = chk_oh2(self.gui.chk_oh2.checkState())
            if chk_oh2:
                grupos.append('chk_oh2')
                score += 2      
            #desativante
            chk_arila = lambda x: True if x == 2 else False
            chk_arila = chk_arila(self.gui.chk_arila.checkState())
            if chk_arila:
                grupos.append('chk_arila')
                score -= 1
            chk_metil = lambda x: True if x == 2 else False
            chk_metil = chk_metil(self.gui.chk_metil.checkState())
            if chk_metil:
                grupos.append('chk_metil')
                score -= 1
            return score, grupos
        except Exception as e:
            self.show_error_message(e)    
       
    def fCalcularPotencia(self, smiles):
        try:
            smiles = self.fLerSmiles()
            valorH1, valorH2, scoreH, grupo_etil_ligado_ao_nitrogenio = self.calc.alfaH_score(smiles)
            if valorH1 or valorH2:
                if grupo_etil_ligado_ao_nitrogenio:
                    self.gui.ln_grupoetil.setText('Grupo Etil ligado ao Nitrogênio')
                self.gui.h_alfa1.setText(str(valorH1))
                self.gui.h_alfa2.setText(str(valorH2))
                self.gui.score_halfalados.setText(str(scoreH))

                scoreA, grupos = self.fLerCheckbox()

                if 'chk_acido' in grupos:
                    self.gui.ln_acidocarb.setText('+3')
                if 'chk_pirrolidina' in grupos:
                    self.gui.ln_pirrolidina.setText('+3')
                if 'chk_anelS' in grupos:
                    self.gui.ln_anelS.setText('+3')
                if 'chk_5ou6' in grupos:
                    self.gui.ln_5ou6.setText('+2')
                if 'chk_morfolina' in grupos:
                    self.gui.ln_morfolna.setText('+1')
                if 'chk_7membros' in grupos:
                    self.gui.ln_7membros.setText('+1')
                if 'chk_5C' in grupos:
                    self.gui.ln_5C.setText('+1')
                if 'chk_retirador1' in grupos:
                    self.gui.score_ret.setText('+1')
                if 'chk_retirador2' in grupos:
                    self.gui.score_ret.setText('+2')
                if 'chk_oh1' in grupos:
                    self.gui.ln_hidroxila.setText('+1')
                if 'chk_oh2' in grupos:
                    self.gui.ln_hidroxila.setText('+2')
                if 'chk_oh2' in grupos:
                    self.gui.ln_hidroxila.setText('+2')
                if 'chk_arila' in grupos:
                    self.gui.ln_arila.setText('-1')
                if 'chk_metil' in grupos:
                    self.gui.ln_metil.setText('-1')      

                potScore = scoreH + scoreA
                if potScore >= 4:
                    cpc = 4
                    ai = 1500
                elif potScore == 3:
                    cpc = 3
                    ai = 400
                elif potScore == 2:
                    cpc = 2
                    ai = 100
                else:
                    cpc = 1
                    ai = 18
                self.gui.score_final.setText(str(potScore))
                self.gui.categoria.setText(str(cpc))
                self.gui.limite_diario.setText(str(ai))  
        except Exception as e:
            self.show_error_message(e)        

## FECHAR PROGRAMA
    def closeEvent(self, event):
        reply = QMessageBox.question(self, 'Fechar programa', 'Você realmente deseja fechar essa aplicação?',
                                    QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            event.accept()        
        else:
            event.ignore()

    def show_error_message(self, exception):
        if isinstance(exception, Exception):  # Verifica se a variável 'exception' é uma instância de Exception
            error_message = f"Erro: {exception}"
        else:
            error_message = "Não foi possível ler o texto."
        QMessageBox.information(None, "Erro de Leitura", error_message, QMessageBox.Ok)


if __name__ == '__main__':
    import sys
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    
    gui = App1()
    gui.show()
    sys.exit(app.exec())
        