from get_knowledge import get_id, get_activity, get_bind, get_text, get_drug
from create_prompt import promptt
from deepseek import chat
import re
import sys
from rdkit import Chem

def extract_smiles(text):
    pattern = r'<smiles>(.*?)</smiles>'
    smiles_list = re.findall(pattern, text)
    if len(smiles_list) > 0:
        return smiles_list[0]
    else:
        return 'not provided'

def is_valid_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0
        # Optional: Kekulization verification (check aromaticity, etc.)
        try:
            Chem.SanitizeMol(mol)
            return 1
        except:
            return 0
    except:
        return 0

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py pref_name1 [pref_name2 ...]")
        sys.exit(1)
    
    pref_name_list = sys.argv[1:]  # Get all command-line arguments after the script name
    
    for i in pref_name_list:
        pref_name = i
        try:
            uniprotid, chemblid = get_id(pref_name)
            activity = get_activity(chemblid)
            drug = get_drug(chemblid)
            binding_site = get_bind(uniprotid)
            text = get_text(pref_name)
            prompt = promptt(pref_name=pref_name, Ki=activity, binding_site=binding_site, 
                            request='1', chemblid=chemblid, text=text)
            prompt_list = prompt.get_prompt()
            cont = 1
            print(chemblid)
            
            while(cont):
                answer = chat(prompt_list)
                smiles = extract_smiles(answer)
                if is_valid_smiles(smiles):
                    print(smiles)
                    cont = 0
                    break
            
            with open('output.txt', 'a') as file:
                file.write(chemblid)
                file.write(' ')
                file.write(smiles)
                file.write('\n')
            print('=====================')
        except Exception as e:
            print(f"Error processing {pref_name}: {str(e)}")
            continue

if __name__ == "__main__":
    main()
