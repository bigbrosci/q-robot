import os
import subprocess

def force_recalc_libreoffice(file_path):
    """
    Open an Excel file in LibreOffice Calc headless mode to force a full recalculation of formulas,
    then save the file. This function uses LibreOffice's --headless and --convert-to options.
    
    Parameters:
      file_path (str): Path to the Excel file.
    
    Returns:
      None
    """
    directory = os.path.dirname(file_path)
    # The "--calc" flag ensures LibreOffice opens the file as a Calc file.
    # "--convert-to xlsx" will save the file in xlsx format, which triggers recalculation.
    command = [
        "soffice",
        "--headless",
        "--calc",
        "--convert-to",
        "xlsx",
        file_path,
        "--outdir",
        directory
    ]
    subprocess.run(command, check=True)
    print(f"Recalculation and update complete for {file_path}")

# Example usage:
# Suppose you have a folder with many Excel files
def process_all_excel_files(root_dir):
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file.endswith('.xlsx'):
                file_path = os.path.join(root, file)
                try:
                    force_recalc_libreoffice(file_path)
                except subprocess.CalledProcessError as e:
                    print(f"Error processing {file_path}: {e}")

# Usage:
process_all_excel_files(".")

