from get_knowledge import get_id,get_activity,get_bind,get_text,get_drug
from create_prompt import promptt
from deepseek import chat
import re
from rdkit import Chem
def extract_smiles(text):
    pattern = r'<smiles>(.*?)</smiles>'
    smiles_list = re.findall(pattern, text)
    if len(smiles_list)>0 :
        return smiles_list[0]
    else:
        return 'not provided'

def is_valid_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0
        # 可选：进行Kekulization验证（检查芳香性等）
        try:
            Chem.SanitizeMol(mol)
            return 1
        except:
            return 0
    except:
        return 0
# pref_name_list = ['Epidermal growth factor receptor erbB1','Coagulation factor VII'] 在此键入你的pref_name即可
pref_name_list = ['Epidermal growth factor receptor erbB1','Coagulation factor VII']
for i in pref_name_list:
    pref_name = i
    uniprotid,chemblid = get_id(pref_name)
    activity = get_activity(chemblid)
    drug = get_drug(chemblid)
    binding_site = get_bind(uniprotid)
    text = get_text(pref_name)
    prompt = promptt(pref_name=pref_name, Ki = activity, binding_site=binding_site, request='1',chemblid=chemblid,text=text)
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
