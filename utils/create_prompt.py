class promptt:
    def __init__(self,pref_name,binding_site,Ki,request,chemblid,text):
        self.pref_name = pref_name
        self.binding_site = binding_site
        self.Ki = Ki
        self.request = request
        self.chemblid = chemblid
        self.text = text

    def rule(self):
        prompt = '''Define the rule of generating SMILES of ligand.
        '''
        return prompt
        
    def get_binding_site_propmt(self):
        prompt = f'''###Background: You are an expert in the field of molecular generation, proficient in generating ligand molecules for proteins based on certain information. Now, you have accepted the literature and summarized some knowledge from it. You can obtain the knowledge you have learned through context. Your current stage: binding site analysis (structural docking), studying the binding site information of the target I provided to you.
###Role: Molecular Generation Expert (Structural Analysis Mode)
###Specific task: At this stage, your task is to analyze some binding site related information of the target downloaded from the uniprotoKB database in JSON format, learn the structure related knowledge of the target from you, and use it as the context for subsequent molecular generation.
###Input: {self.binding_site} (UniProt exported binding site data)
###Analysis steps:
1. Extract structural features that combine positions
2. Identify the molecules or residues that can bind to the binding site and analyze their characteristics
3. Try to piece together the molecular residues bound to that position into a single molecule
4. Macroscopically provide the spatial significance of the binding site of the target, laying the foundation for the subsequent generation of bound molecules.
(The above steps are guesses from non chemistry students. If you have better analytical steps, please make your own decision, in order to better learn the knowledge of the binding site of the target mentioned above.)
###Output: Knowledge learned
'''
        return prompt


    def get_resorted_prompt(self):
        # 将文本内容放在json材料后
        prompt1 = self.get_article()
        prompt2 = self.get_binding_site_propmt()
        prompt3 = self.get_Ki_ligand()
        prompt4 = self.get_generate_prompt()
        prompt5 = self.get_check_prompt()
        rule = self.rule()
        return rule,prompt2,prompt3,prompt1,prompt4,prompt5,self.chemblid

    def get_icnm(self,icnm):
        prompt = f'''###Background: You are an expert in the field of molecular generation, proficient in generating ligand molecules for proteins based on certain information. You can obtain the knowledge you have learned through context. Your current stage: Experimental validation has been combined with ligand analysis to study the information on the existing binding ligands of the target that I have provided to you.
###Role: Molecular Generation Expert (Structural Analysis Mode)
###Specific task: At this stage, your task is to analyze some highly active ligand related information of the target downloaded from the CHEMBL database in JSON format, learn about the ligand related knowledge of the target from you, and use it as the context for subsequent molecule generation.
###Input: {icnm} (icnm binding site data exported from UniProt)
###Analysis steps:
1. Establish a 2D-QSAR model: identify key descriptors that affect pIC50 (such as AlogP, TPSA)
2. Conduct 3D pharmacophore comparison: Calculate molecular overlap RMSD and extract conservative features
3. Building an active cliff analysis: identifying substituent patterns that lead to sudden changes in activity
(The above steps are guesses from non chemistry students. If you have better analytical steps, please make your own decision, with the aim of better learning about the relevant ligands of the target mentioned above.)
###Output: Knowledge learned
'''
        return prompt

    def get_icnm_prompt(self,icnm):
        prompt1 = self.get_article()
        prompt2 = self.get_binding_site_propmt()
        prompt3 = self.get_icnm(icnm)
        prompt4 = self.get_generate_prompt()
        prompt5 = self.get_check_prompt()
        rule = self.rule()
        return rule,prompt1,prompt2,prompt3,prompt4,prompt5,self.chemblid

    def get_article(self):
        prompt = f'''###Background: You are an expert in the field of molecular generation and can proficiently generate ligand molecules based on some information of proteins. Smiles said that you will gain relevant knowledge step by step through multiple rounds of conversations with me, and gradually improve your ability to generate target ligand molecules. Your current stage: literature knowledge structuring (knowledge anchoring), you can conduct ligand related research on a certain protein target based on the following literature information about the target.
###Role: Molecular Generation Expert (Literature Analysis Mode)
###Input: {self.text} (target information from paper/network)
###Task: Perform the following analysis process:
1. Basic information analysis of the target
2. Analysis of structural information of the target
3. Analysis of experimental information on the known target
4. Analysis of ligand generation mode of the target
(The above steps are guesses from non chemistry students. If you have better analytical steps, please make your own decision, with the aim of better learning the knowledge of the target article mentioned above.)
###Output: Knowledge obtained from the article.
        '''
        return prompt
    
    def get_Ki_ligand(self):
        prompt = f'''
###Background: You are an expert in the field of molecular generation, proficient in generating ligand molecules for proteins based on certain information. You can generate the smiles representation of the ligand molecule for a certain protein target based on the following literature information, binding sites, existing drug molecules, and active substances with high binding affinity. Your current stage: ligand SAR analysis (structure-activity relationship modeling)
###Task: At this stage, your task is to analyze the ligand related information of some high binding affinity values (arranged in descending order of binding ability) of the target downloaded from the Chembl database in JSON format, and learn the relevant knowledge of the target from you as the context for subsequent molecule generation.
###Role: Molecular Generation Expert (QSAR Analysis Mode)
###Input: {self.Ki} (ChEMBL High Activity Molecular Dataset)
###Analysis process:
1. Extract common features of ligands
2. Analyze the binding mode between ligands and targets
3. Extract the molecular fragments that make up the ligand in the example
(The above steps are guesses from non chemistry students. If you have better analytical steps, please make your own decision, with the aim of better learning about the relevant ligands of the target mentioned above.)
###Output: Knowledge learned
        '''
        return prompt
    
    def get_generate_prompt(self):
        prompt = f'''###Background: You are an expert in the field of molecular generation, proficient in generating ligand molecules for proteins based on certain information. You can generate the smiles representation of the ligand molecule for a certain protein target based on the following literature information, binding sites, existing drug molecules, and active substances with high binding affinity. Current stage: Molecular generation (rational design), designing ligand molecules for the target.
###Task: At this stage, your task is to generate new ligand molecules for the target based on the knowledge you have obtained, in order to obtain molecules with higher binding affinity, higher drug like properties, and higher molecular synthesizability compared to the provided knowledge. Please generate the ligand molecule for the target based on the knowledge of "characteristics of high affinity molecules" and "methods for generating high affinity molecules smiles" in the context, think step by step, provide the reasoning process, check the reasoning process, and summarize the final smile results. Require an affinity Vina score below -9.
###Role: Molecular Generation Expert (de novo design pattern)
###Molecular design logic:
<Step 1>Find the core structure in the ligand as the reference molecule.
<Step 2>Introduce drug-resistant groups (such as rigid macrocyclic structures).
<Step 3>Optimize the molecule for easier synthesis.
Ensure the generation of legal and high binding molecules.
###Output specification:
xml
<design-process>
<base>Reference molecule</base>
<step1>
<modify>First modified content</modify>
<smiles>First modified smiles</smiles>
</step1>
<step2>
<modify>Second modification content</modify>
<smiles>Second modified smiles</smiles>
</step2>
<final>Final smiles</final>
</design-process>
###Attention: Max Ring Size is smaller than 8.
        '''
        return prompt
    
    def get_check_prompt(self):
        prompt = f'''Background: You are an expert in the field of molecular generation and can proficiently generate ligand molecules for proteins based on certain information. You can generate the smiles representation of the ligand molecule for a certain protein target based on the following literature information, binding sites, existing drug molecules, and active substances with high binding affinity. Your current stage is molecular validation (multidimensional evaluation), where you need to evaluate the molecules generated in the previous section, validate, rank, and output them.
###Task: At this stage, your task is to carefully read the molecules you generate in the context, perform verification, validation, and ranking work, and ensure that the smiles are valid. It is forbidden to generate invalid molecules!!! Definitely generate molecules with strong binding affinity!!!! Among them, the molecule ranked first must be the one with the strongest binding force, and we must spare no effort to improve its binding ability, while maintaining its legitimacy.
###Verification Protocol:
1. Chemical validity check: Use RDKit to verify the correctness of SMILES syntax/valence state.
2. Docking verification: AutoDock Vina rapid docking with the original crystal structure (Δ G ≤ -9 kcal/mol)
3. Drug Evaluation: QED ≥ 0.6 and SAscore ≤ 4
4. Comprehensive sorting: Weighted scores based on (0.4 × pKiupred)+(0.3 × QED)+(0.3 × SAscore)
###Output template:
xml
<validation-summary>
<smiles>The only legitimate and highly binding molecule for testing</smiles>
<qed>Predict its qed indicators</qed>
<affinity>Predict their Vina score</affinity>
<sa>Predict its synthesizability index sa</sa>
</validation-summary>
###Attention: The entire process must ensure the generation of legal smiles!!!!!!!!!! The conditions for legal smiles:
1. Grammar level: Atomic symbols must be valid (such as C, [Na+]), bond connections must be clear (single bonds can be omitted, double/triple bonds use=#), parentheses are closed in pairs, ring markers match numbers (such as C1CCCC1), aromatic atoms are lowercase (such as benzene c1ccccc1); 2. At the chemical level: the atomic valence state is reasonable (such as carbon not exceeding 4 bonds), the hydrogen atom is implicitly expressed (special explicit hydrogen such as [NH3]), the charge labeling is correct (such as [Fe+2]), and the overall structure is free of contradictions (such as ring energy closure and ion charge balance).
###Attention: Please try to avoid the generation of molecules in the following situations: ① Large rings with more than 8 molecules appear The situation where three rings share a molecule ③ There are too many heteroatoms in the molecule. Heteroatoms in the molecule appear multiple times in the same ring.
'''
        return prompt
    
    def get_prompt(self):
        prompt1 = self.get_article()
        prompt2 = self.get_binding_site_propmt()
        prompt3 = self.get_Ki_ligand()
        prompt4 = self.get_generate_prompt()
        prompt5 = self.get_check_prompt()
        # smiles = self.rule()
        return prompt1,prompt2,prompt3,prompt4,prompt5,self.chemblid

    def get_prompt_from_id(self,uniprotid,chemblid):
        # 设置 URL
        url = f"https://www.ebi.ac.uk/chembl/api/data/target/{chemblid}.json"
        
        # 发送 GET 请求
        response = requests.get(url)
        
        # 检查请求是否成功
        if response.status_code == 200:
            # 解析 JSON 数据
            data = response.json()
            return data['pref_name']
        else:
            return None

    def get_function(self,function,paper):
        prompt = f'''Hello, you are an experienced ligand generation expert for protein targets. Now I will give you the functional information of the protein target {self.pref_name}, please study this information: {function}. After studying the structural information, please read the current relevant literature: {paper}. Alright, after reading the literature and experiencing the functional information of this target, please try generating the molecular smiles representation of the ligand that binds to this target. The final result will be wrapped in<smiles></smiles>for subsequent extraction.
        '''
    def get_no_rag(self):
        prompt = f'''Please generate the ligand molecule smiles for the {self.pref_name} target and wrap them with<smiles></smiles>for subsequent extraction'''
        return prompt

        
