import subprocess
import os

def install_pyqt5():
    powershell_code = r'''
$package = "PyQt5"
$packageInstalled = Get-Module -ListAvailable | Where-Object { $_.Name -eq $package }

if (-not $packageInstalled) {
    Write-Host "PyQt5 não está instalado. Instalando..."
    & pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org PyQt5 rdkit pandas
} else {
    Write-Host "PyQt5 já está instalado."
}
    '''

    # Salva o código PowerShell em um arquivo temporário
    ps_script_path = os.path.join(os.environ["TEMP"], "install_pyqt5.ps1")
    with open(ps_script_path, "w") as ps_file:
        ps_file.write(powershell_code)

    try:
        # Executa o arquivo PowerShell com powershell.exe
        subprocess.run(["powershell.exe", "-ExecutionPolicy", "Bypass", "-File", ps_script_path], shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Erro ao executar o PowerShell: {e}")
    else:
        print("PyQt5 instalado com sucesso!")

if __name__ == "__main__":
    install_pyqt5()
