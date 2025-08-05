import requests
import json
import time
def download_alphafold_pdb(chembl_id):
    # Step 1: 从ChEMBL获取UniProt ID
    chembl_url = f"https://www.ebi.ac.uk/chembl/api/data/target/{chembl_id}?format=json"
    response = requests.get(chembl_url)
    if response.status_code != 200:
        print(f"错误：未找到ChEMBL ID {chembl_id}")
        return

    data = response.json()
    uniprot_id = data.get('target_components', [{}])[0].get('accession')
    if not uniprot_id:
        print("错误：未找到关联的UniProt ID")
        return

    # Step 2: 下载AlphaFold PDB
    af_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    pdb_response = requests.get(af_url)
    if pdb_response.status_code == 200:
        with open(f"/root/autodl-tmp/pdb_files/{chembl_id}.pdb", "wb") as f:
            f.write(pdb_response.content)
        print(f"文件已保存为 AF-{uniprot_id}.pdb")
    else:
        print(f"错误：{chemblid}AlphaFold模型不存在")

if '__name__' == '__main__':
    with open('full_crossdock_deepseek_prompt.jsonl') as f:
        all_data = [json.loads(line) for line in f]    
    intt = 0
    for i in all_data[900:]:
        intt = intt + 1
        download_alphafold_pdb(i['chemblid'])
        if intt%10 == 0:
            time.sleep(10)
    # for i in data:
    #     download_alphafold_pdb(i[-1])
