import numpy as np
from typing import List, Optional, Tuple

def evidence_theoretic_score(
    v_d: float,
    qed: float,
    sa: float,
    constraints: Optional[List[float]] = None,
    tau_g: float = -8.0,
    tau_b: float = -6.0,
    k_g: float = 0.5,
    k_b: float = 0.5,
    alpha: float = 1.5
) -> Tuple[float, dict]:
    """
    Compute evidence-theoretic score for ligand evaluation.
    
    Args:
        v_d: Vina docking score (more negative = stronger binding)
        qed: QED score (0-1, higher is better)
        sa: SA score (0-1, higher is better)
        constraints: List of constraint satisfaction factors (0-1, optional)
        tau_g: Threshold for strong binding (supporting 'good' hypothesis)
        tau_b: Threshold for weak binding (supporting 'bad' hypothesis)
        k_g: Steepness parameter for good hypothesis sigmoid
        k_b: Steepness parameter for bad hypothesis sigmoid
        alpha: Exponent for QED power transformation
    
    Returns:
        Tuple containing:
            - Final score (0-100)
            - Dictionary with belief, plausibility, and uncertainty
    """
    
    # Helper function: sigmoid for evidence assignment
    def sigmoid(x, tau, k):
        return 1.0 / (1.0 + np.exp(-k * (x - tau)))
    
    # Helper function: combine two BPAs using Dempster's rule
    def combine_bpas(bpa1, bpa2):
        """Combine two Basic Probability Assignments."""
        K = (bpa1['G'] * bpa2['B']) + (bpa1['B'] * bpa2['G'])
        
        if K >= 1.0:
            return {'G': 0.0, 'B': 0.0, 'GB': 1.0}
        
        m_g = (bpa1['G'] * bpa2['G'] + bpa1['G'] * bpa2['GB'] + bpa1['GB'] * bpa2['G']) / (1.0 - K)
        m_b = (bpa1['B'] * bpa2['B'] + bpa1['B'] * bpa2['GB'] + bpa1['GB'] * bpa2['B']) / (1.0 - K)
        m_gb = (bpa1['GB'] * bpa2['GB']) / (1.0 - K)
        
        # Normalize
        total = m_g + m_b + m_gb
        return {'G': m_g/total, 'B': m_b/total, 'GB': m_gb/total}
    
    # 1. Compute Basic Probability Assignments for each metric
    # Docking score evidence
    m_dock_g = sigmoid(v_d, tau_g, k_g)
    m_dock_b = sigmoid(-v_d, tau_b, k_b)
    if m_dock_g + m_dock_b > 1.0:
        m_dock_g /= (m_dock_g + m_dock_b)
        m_dock_b /= (m_dock_g + m_dock_b)
        m_dock_gb = 0.0
    else:
        m_dock_gb = 1.0 - m_dock_g - m_dock_b
    
    bpa_dock = {'G': m_dock_g, 'B': m_dock_b, 'GB': m_dock_gb}
    
    # QED evidence
    m_qed_g = qed ** alpha
    m_qed_b = (1.0 - qed) ** alpha
    if m_qed_g + m_qed_b > 1.0:
        m_qed_g /= (m_qed_g + m_qed_b)
        m_qed_b /= (m_qed_g + m_qed_b)
        m_qed_gb = 0.0
    else:
        m_qed_gb = 1.0 - m_qed_g - m_qed_b
    
    bpa_qed = {'G': m_qed_g, 'B': m_qed_b, 'GB': m_qed_gb}
    
    # SA evidence (deterministic)
    bpa_sa = {'G': sa, 'B': 1.0 - sa, 'GB': 0.0}
    
    # 2. Combine BPAs sequentially
    bpa_combined = combine_bpas(bpa_dock, bpa_qed)
    bpa_combined = combine_bpas(bpa_combined, bpa_sa)
    
    # 3. Add constraints if provided
    if constraints:
        for constraint_satisfaction in constraints:
            bpa_constraint = {'G': constraint_satisfaction, 'B': 1.0-constraint_satisfaction, 'GB': 0.0}
            bpa_combined = combine_bpas(bpa_combined, bpa_constraint)
    
    # 4. Compute belief, plausibility, and final score
    belief = bpa_combined['G']
    plausibility = bpa_combined['G'] + bpa_combined['GB']
    final_score = 100.0 * (belief + plausibility) / 2.0
    
    # 5. Return result
    return final_score, {
        'belief': belief,
        'plausibility': plausibility,
        'uncertainty': bpa_combined['GB'],
        'confidence_interval': [belief, plausibility]
    }


# 使用示例
if __name__ == "__main__":
    # 示例1: 优秀配体
    score1, info1 = evidence_theoretic_score(v_d=-9.5, qed=0.85, sa=0.90)
    print(f"优秀配体: {score1:.1f}/100, 不确定性: {info1['uncertainty']:.3f}")
    
    # 示例2: 冲突证据
    score2, info2 = evidence_theoretic_score(v_d=-8.2, qed=0.30, sa=0.25)
    print(f"冲突配体: {score2:.1f}/100, 不确定性: {info2['uncertainty']:.3f}")
    
    # 示例3: 带有约束的配体
    score3, info3 = evidence_theoretic_score(
        v_d=-7.0, qed=0.65, sa=0.60, 
        constraints=[0.8, 1.0]
    )
    print(f"约束配体: {score3:.1f}/100, 置信区间: {info3['confidence_interval']}")
    
    # 批量处理示例
    ligands = [
        {"v_d": -8.5, "qed": 0.75, "sa": 0.80},
        {"v_d": -6.5, "qed": 0.90, "sa": 0.85},
        {"v_d": -9.0, "qed": 0.45, "sa": 0.50},
    ]
    
    print("\n批量处理结果:")
    for i, ligand in enumerate(ligands, 1):
        score, info = evidence_theoretic_score(**ligand)
        print(f"配体{i}: {score:.1f}, 不确定性: {info['uncertainty']:.3f}")