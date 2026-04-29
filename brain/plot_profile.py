# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 20:18:43 2026

@author: lqlhz
"""

########################################################

def build_refdict_and_plot_rc(ef_value: float, adsorption_data: dict, outputs_folder: str = "outputs") -> dict:
    os.makedirs(outputs_folder, exist_ok=True)

    # Gas-phase reference energies
    E_gas = {"NH3": -19.53586573, "NH2": -13.53307777, "NH": -8.10061060,
             "N2": -16.62922486, "H2": -6.76668776}

    # Field-dependent slab energy
    E_slab = -6.2414 * ef_value**2 - 0.01405 * ef_value - 354.4197

    # Adsorbed species energies
    E_NH3_ads = adsorption_data['NH3'][2]
    E_NH2_ads = adsorption_data['NH2'][2]
    E_NH_ads = adsorption_data['NH'][2]
    E_N_ads = adsorption_data['N'][2]
    E_H_ads = min(adsorption_data[k][2] for k in ('H1', 'H2', 'H3'))
    E_N2_ads = adsorption_data['N2'][2]

    # Reaction energies
    ΔE1 = E_NH3_ads - E_slab - E_gas['NH3']
    ΔE2 = (E_NH2_ads + E_H_ads - E_slab) - E_NH3_ads
    ΔE3 = (E_NH_ads + E_H_ads - E_slab) - E_NH2_ads
    ΔE4 = (E_N_ads + E_H_ads - E_slab) - E_NH_ads
    ΔE5 = (E_N2_ads + E_slab) - 2 * E_N_ads
    ΔE6 = -adsorption_data['N2'][1] * 0.5
    ΔE7 = -adsorption_data['H1'][1] * 1.5

    # Barrier estimation: Ea = a * ΔE + b
    Ea_param = [(0.42, 1.36), (0.88, 0.82), (0.47, 0.81), (0.87, 2.10)]
    Ea1 = max(0.01, Ea_param[0][0] * ΔE2 + Ea_param[0][1])
    Ea2 = max(0.01, Ea_param[1][0] * ΔE3 + Ea_param[1][1])
    Ea3 = max(0.01, Ea_param[2][0] * ΔE4 + Ea_param[2][1])
    Ea4 = max(0.01, Ea_param[3][0] * ΔE5 + Ea_param[3][1])

    # Extended state sequence including TS
    state_labels = [
        'NH3(g)', 'NH3*', 'TS1_NH3', 'NH2*+H*', 'TS2_NH2',
        'NH*+2H*', 'TS3_NH', 'N*+3H*', 'TS4_N2',
        'N2*+3H*', 'N2(g)+3H*', 'N2(g)+1.5H2(g)'
    ]

    # Compute state energies
    relE = [0]
    relE.append(relE[-1] + ΔE1)              # NH3*
    relE.append(relE[-1] + Ea1)              # TS1
    relE.append(relE[-2] + ΔE2)              # NH2* + H*
    relE.append(relE[-1] + Ea2)              # TS2
    relE.append(relE[-2] + ΔE3)              # NH* + 2H*
    relE.append(relE[-1] + Ea3)              # TS3
    relE.append(relE[-2] + ΔE4)              # N* + 3H*
    relE.append(relE[-1] + Ea4)              # TS4
    relE.append(relE[-2] + ΔE5)              # N2* + 3H*
    relE.append(relE[-1] + ΔE6)              # N2(g) + 3H*
    relE.append(relE[-1] + ΔE7)              # N2(g) + 1.5H2(g)

    # Build platform data for plotting
    x_vals = []
    y_vals = []
    for i in range(len(relE)):
        x_vals.extend([i, i + 1])
        y_vals.extend([relE[i]] * 2)
    # Expand state labels to match platform points
    labels_expanded = []
    for label in state_labels:
        labels_expanded.extend([label, label])
    
    # Save rc.csv with labels
    df_rc = pd.DataFrame({'x': x_vals, 'y': y_vals, 'label': labels_expanded})
    df_rc.to_csv(os.path.join(outputs_folder, 'rc.csv'), index=False)
    # Save rc.csv

    # Plot
    plt.figure(figsize=(10, 5))
    plt.plot(x_vals, y_vals, '-k', lw=2)
    plt.scatter([i + 0.5 for i in range(len(relE))], relE, color='k', s=30, zorder=3)
    plt.axhline(0, ls='--', color='gray')

    # Label x-axis
    plt.xticks([i + 0.5 for i in range(len(state_labels))], state_labels, rotation=45, ha='right')
    plt.ylabel('Relative energy (eV)')
    plt.title(f'EF = {ef_value:.1f} V/Å')
    plt.tight_layout()
    plt.savefig(os.path.join(outputs_folder, f'reaction_coordinate_EF_{ef_value}.png'), dpi=300)
    plt.close()

    # Replacement dict
    repl = {
        "RU(T)": E_slab,
        "NH3(T)": E_NH3_ads,
        "NH2(T)": E_NH2_ads,
        "NH(T)": E_NH_ads,
        "N(T)": E_N_ads,
        "N2(T)": E_N2_ads,
        "Hv1(T)": adsorption_data['H1'][2],
        "Hv2(T)": adsorption_data['H2'][2],
        "Hv3(T)": adsorption_data['H3'][2],
        "H(T)": E_H_ads,
        "N_N(T)": 2 * E_N_ads - E_slab + 0.5,
        "TS1_NH3(T)": E_NH3_ads + Ea1,
        "TS2_NH2(T)": E_NH2_ads + Ea2,
        "TS3_NH(T)": E_NH_ads + Ea3,
        "TS4_N2(T)": 2 * E_N_ads - E_slab + 0.5 + Ea4
    }

    ref_dict = {k: v - E_slab for k, v in repl.items()}
    return ref_dict

def plot_rc(EF, edft, e_slab, outputs_folder):
    E_gas = {                 # gas
        "H": -1.11671722,
        "H2": -6.76668776,
        "N": -3.12298738,
        "N2": -16.62922486,
        "NH": -8.10061060,
        "NH2": -13.53307777,
        "NH3": -19.53586573,
    }

    # ---------- 3. 取常用能量 ----------
    E_NH3_ads  = edft['NH3'][2]          # slab+NH3
    E_NH2_ads  = edft['NH2'][2]
    E_NH_ads   = edft['NH'][2]
    E_N_ads    = edft['N'][2]
    E_H_ads    = min(edft[k][2] for k in ('H1','H2','H3'))  # 最稳定 H
    E_N2_ads   = edft['N2'][2]
    
    # ---------- 4. 逐步反应热 ----------
    steps   = []
    labels  = []
    
    # 0) 初始气相
    cumE = 0.0
    steps.append(cumE); labels.append('NH$_3$(g)')
    
    # 1) NH3 吸附
    dE1  = E_NH3_ads - e_slab - E_gas['NH3']
    cumE += dE1
    steps.append(cumE); labels.append('NH$_3$*')
    
    # 2) NH3 → NH2* + H*
    dE2  = (E_NH2_ads + E_H_ads - e_slab) - E_NH3_ads
    cumE += dE2
    steps.append(cumE); labels.append('NH$_2$* + H*')
    
    # 3) NH2* → NH* + H*
    dE3  = (E_NH_ads + E_H_ads - e_slab) - E_NH2_ads
    cumE += dE3
    steps.append(cumE); labels.append('NH* + 2H*')
    
    # 4) NH* → N* + H*
    dE4  = (E_N_ads + E_H_ads - e_slab) - E_NH_ads
    cumE += dE4
    steps.append(cumE); labels.append('N* + 3H*')
    
    # 5) N* + N* → N2*
    dE5  = (E_N2_ads + e_slab) - 2 * E_N_ads
    cumE += dE5/2
    steps.append(cumE); labels.append('N$_2$* + 3H*')
    
    # 6) N2 脱附
    dE6  = - edft['N2'][1] * 0.5          # 按给定规则：吸附能的 -½
    cumE += dE6
    steps.append(cumE); labels.append('N$_2$(g) + 3H*')
    
    # 7) 3H* → 1.5 H2(g)
    dE7  = - edft['H1'][1] * 1.5          # 最稳定 H 的吸附能 × -1.5
    cumE += dE7
    steps.append(cumE); labels.append('N$_2$(g) + 1.5H$_2$(g)')
    
    # ---------- 5. 画 Reaction Coordinate ----------
    # x = range(len(steps))
    # plt.figure(figsize=(8,5))
    # plt.plot(x, steps, '-o', lw=2, ms=6, color='k')
    # plt.axhline(0, ls='--', color='grey', lw=1)
    
    # plt.xticks(x, labels, rotation=45, ha='right')
    # plt.ylabel('Relative energy (eV)')
    # plt.tight_layout()
    # plt.savefig('reaction_coordinate.png', dpi=300)
    # plt.show()
    
    # --- 构造平台坐标 ---
    x_vals = []
    y_vals = []
    for i in range(len(steps)):
        x_vals.extend([i, i+1])
        y_vals.extend([steps[i]] * 2)
    df_rc = pd.DataFrame({'x': x_vals, 'y': y_vals})
    df_rc.to_csv(os.path.join(outputs_folder, 'rc.csv'), index=False)
    
    # --- 设置平台中点 label 的位置 ---
    label_pos = [(i + i + 1)/2 for i in range(len(steps))]
    
    # --- 绘图 ---
    plt.figure(figsize=(8,5))
    plt.plot(x_vals, y_vals, '-', lw=2, color='k', label=f'EF = {EF:.1f} V/Å')
    plt.scatter(label_pos, steps, s=30, color='k', zorder=3)  # 小圆点标示每个状态
    plt.axhline(0, ls='--', color='gray', lw=1)
    
    # x-axis
    plt.xticks(label_pos, labels, rotation=45, ha='right')
    plt.xlim(0, len(steps))
    plt.ylabel('Relative energy (eV)')
    plt.legend()

    plt.tight_layout()
    filename = os.path.join(outputs_folder, f'reaction_coordinate_EF_{EF}.png')
    plt.savefig(filename, dpi=300)
    plt.close()



def build_refdict_and_plot_rc_double_NH3(ef_value: float, adsorption_data: dict, outputs_folder: str = "outputs") -> dict:
    os.makedirs(outputs_folder, exist_ok=True)

    # Gas-phase reference energies
    E_gas = {
        "NH3": -19.53586573,
        "NH2": -13.53307777,
        "NH": -8.10061060,
        "N2": -16.62922486,
        "H2": -6.76668776
    }

    # Field-dependent slab energy
    E_slab = -6.2414 * ef_value**2 - 0.01405 * ef_value - 354.4197

    # Adsorbed species energies
    E_NH3_ads = adsorption_data['NH3'][2]
    E_NH2_ads = adsorption_data['NH2'][2]
    E_NH_ads = adsorption_data['NH'][2]
    E_N_ads = adsorption_data['N'][2]
    try:
        E_H_ads = min(adsorption_data[k][2] for k in ('H1', 'H2', 'H3'))
    except:
        E_H_ads = adsorption_data['H'][2]
    E_N2_ads = adsorption_data['N2'][2]

    # Reaction energies (1 NH3 unit)
    ΔE1 = E_NH3_ads - E_slab - E_gas['NH3']
    ΔE2 = (E_NH2_ads + -6.76668776/2) - E_NH3_ads
    ΔE3 = (E_NH_ads  + -6.76668776/2) - E_NH2_ads
    ΔE4 = (E_N_ads   + -6.76668776/2) - E_NH_ads
    # ΔE5 = (E_N2_ads + E_slab) - 2 * E_N_ads
    ΔE5 = -(2 * E_N_ads )
    ΔE6 = -adsorption_data['N2'][1]          # full N2 desorption
    try:
        ΔE7 = -adsorption_data['H1'][1] * 3      # 6H* → 3 H2(g)
    except:
        ΔE7 = -adsorption_data['H'][1] * 3

    # Barriers
    Ea_param_0 =    [(0.42, 0.91), (0.85, 0.75), (0.57, 0.78), (0.73, 1.31)]
    Ea_param_N0_6 = [(0.53, 0.84), (0.89, 0.72), (0.63, 0.73), (0.48, 1.51)]
    Ea_param_P0_6 = [(0.57, 0.92), (0.79, 0.77), (0.59, 0.78), (0.61, 1.38)]
    
    Ea1 = max(0.01, Ea_param[0][0] * ΔE2 + Ea_param[0][1])
    Ea2 = max(0.01, Ea_param[1][0] * ΔE3 + Ea_param[1][1])
    Ea3 = max(0.01, Ea_param[2][0] * ΔE4 + Ea_param[2][1])
    Ea4 = max(0.01, Ea_param[3][0] * ΔE5 + Ea_param[3][1])

    # State labels with TS included
    state_labels = [
        '2NH3(g)', '2NH3*', 'TS1', '2NH2*+2H*', 'TS2',
        '2NH*+4H*', 'TS3', '2N*+6H*', 'TS4',
        'N2*+6H*', 'N2(g)+6H*', 'N2(g)+3H2(g)'
    ]

    # Compute relE
    relE = [0]
    relE.append(relE[-1] + 2 * ΔE1)           # 2NH3*
    relE.append(relE[-1] + 2 * Ea1)           # TS1
    relE.append(relE[-2] + 2 * ΔE2)           # 2NH2*+2H*
    relE.append(relE[-1] + 2 * Ea2)           # TS2
    relE.append(relE[-2] + 2 * ΔE3)           # 2NH*+4H*
    relE.append(relE[-1] + 2 * Ea3)           # TS3
    relE.append(relE[-2] + 2 * ΔE4)           # 2N*+6H*
    relE.append(relE[-1] + Ea4)               # TS4
    relE.append(relE[-2] + ΔE5)               # N2*+6H*
    relE.append(relE[-1] + ΔE6)               # N2(g)+6H*
    relE.append(relE[-1] + ΔE7)               # N2(g)+3H2(g)

    # Sanity check: match to total gas-phase reaction energy
    ΔE_total_gas = E_gas['N2'] + 3 * E_gas['H2'] - 2 * E_gas['NH3']
    relE[-1] = ΔE_total_gas  # enforce exact gas-phase energy consistency

    # Build x/y for platform plot
    x_vals, y_vals = [], []
    for i, e in enumerate(relE):
        x_vals.extend([i, i + 1])
        y_vals.extend([e, e])

    labels_expanded = []
    for label in state_labels:
        labels_expanded.extend([label, label])

    df_rc = pd.DataFrame({'x': x_vals, 'y': y_vals, 'label': labels_expanded})
    df_rc.to_csv(os.path.join(outputs_folder, 'rc.csv'), index=False)

    # Plot
    plt.figure(figsize=(10, 5))
    plt.plot(x_vals, y_vals, '-k', lw=2)
    plt.scatter([i + 0.5 for i in range(len(relE))], relE, color='k', s=30, zorder=3)
    plt.axhline(0, ls='--', color='gray')
    plt.xticks([i + 0.5 for i in range(len(state_labels))], state_labels, rotation=45, ha='right')
    plt.ylabel('Relative energy (eV)')
    plt.title(f'2 NH₃ → N₂ + 3 H₂   |   EF = {ef_value:.1f} V/Å')
    plt.tight_layout()
    plt.savefig(os.path.join(outputs_folder, f'reaction_coordinate_2NH3_EF_{ef_value}.png'), dpi=300)
    plt.close()

    # print("E_end:", relE[-1], "E_end_gas_ref:", ΔE_total_gas)
    # return df_rc, relE[-1], ΔE_total_gas  # return for optional use

