#############################################################################
# 文件: G2_Extended_Subsystem_Lib.g
# 描述: G2(2) 实形式扩展图分析与子系统筛选 (函数库版)
#############################################################################

# =============================================================================
# 第一部分: 基础数学库 (G2 根系与 Vogan 图)
# =============================================================================

# 1. G2 定义
# 简单根采用 (a1, a2) 坐标表示，其中 a1 为短根方向、a2 为长根方向
SimpleRoots := [[1, 0], [0, 1]];;

# 内积矩阵
# 约定 (a1,a1)=2, (a2,a2)=6, (a1,a2)=(a2,a1)=-3
G2_InnerProdMatrix := [[2, -3], [-3, 6]];;

# 内积函数
# 计算 (v1, v2) = v1^T * M * v2
InnerProduct := function(v1, v2)
    return v1[1] * G2_InnerProdMatrix[1][1] * v2[1] +
           v1[1] * G2_InnerProdMatrix[1][2] * v2[2] +
           v1[2] * G2_InnerProdMatrix[2][1] * v2[1] +
           v1[2] * G2_InnerProdMatrix[2][2] * v2[2];
end;;

# F4(4) Extended Vogan Diagram 子系统生成
# 基于 Weyl 群变换，可按给定位置与涂黑目标筛选配置
FindF4_B4_Subsystems := function(custom_roots, custom_colors)
    local F4_CartanMatrix, F4_InnerProdMatrix, F4_InnerProduct, F4_Reflection,
          F4_IsRootBlack, F4_FormatRoot, F4_ExtendedRoots,
          results, idx, w_roots, k, i, black_count, white_count, colors_list, type_info,
          queue, seen, states, state, key, MakeKey,
          m, F4_IsValidSimpleCoeffRootVec, coeff_bounds,
          all_nonneg, all_nonpos, valid_state,
          all_pos_labels, selected_pos_labels, selected_indices, label, p,
          target_black_labels, required_black_count, target_white_count,
          black_labels_actual, is_match, target_colors, target_color_patterns,
          BuildSelectedIndices, BuildStates, BuildRecord;
    
    # F4 Cartan 矩阵
    F4_CartanMatrix := [
        [ 2, -1,  0,  0],  # a1
        [-1,  2, -2,  0],  # a2
        [ 0, -1,  2, -1],  # a3 (注意: a2-a3 是双键)
        [ 0,  0, -1,  2]   # a4
    ];

    F4_InnerProdMatrix := [
        [4, -2, 0, 0],
        [-2, 4, -2, 0],
        [0, -2, 2, -1],
        [0, 0, -1, 2]
    ];

    F4_InnerProduct := function(v1, v2)
        local a, b, res;
        res := 0;
        for a in [1..4] do
            for b in [1..4] do
                res := res + v1[a] * F4_InnerProdMatrix[a][b] * v2[b];
            od;
        od;
        return res;
    end;
    
    # F4 简单反射（在 simple-root 系数坐标中）
    # s_i(beta) = beta - <beta, alpha_i^vee> * alpha_i
    # 其中 <beta, alpha_i^vee> = sum_j beta_j * A_{j,i}
    F4_Reflection := function(v, alpha_idx)
        local i, result;
        result := ShallowCopy(v);
        m := 0;
        for i in [1..4] do
            m := m + v[i] * F4_CartanMatrix[i][alpha_idx];
        od;
        result[alpha_idx] := result[alpha_idx] - m;
        return result;
    end;
    
    # F4(4) 分裂实形式的根颜色判断
    # 规则: a1 是非紧的 (涂黑)，其他根是紧的
    F4_IsRootBlack := function(root)
        local coeff_a1;
        # 将根表示为简单根的线性组合
        # 如果 a1 的系数为奇数，则为非紧 (涂黑)
        coeff_a1 := root[1];  # 假设 root[1] 是 a1 的系数
        return (coeff_a1 mod 2) <> 0;
    end;
    
    # F4 根格式化
    F4_FormatRoot := function(v)
        local s, i, first;
        s := "";
        first := true;
        for i in [1..4] do
            if v[i] <> 0 then
                if not first and v[i] > 0 then s := Concatenation(s, "+"); fi;
                if v[i] = 1 then s := Concatenation(s, "a", String(i));
                elif v[i] = -1 then s := Concatenation(s, "-a", String(i));
                else s := Concatenation(s, String(v[i]), "a", String(i)); fi;
                first := false;
            fi;
        od;
        if s = "" then s := "0"; fi;
        return s;
    end;
    
    coeff_bounds := [2, 3, 4, 2];
    
    F4_IsValidSimpleCoeffRootVec := function(v)
        local j;
        if not IsList(v) or Length(v) <> 4 then
            return false;
        fi;
        all_nonneg := true;
        all_nonpos := true;
        for j in [1..4] do
            if not IsInt(v[j]) then return false; fi;
            if v[j] > 0 then all_nonpos := false; fi;
            if v[j] < 0 then all_nonneg := false; fi;
            if AbsInt(v[j]) > coeff_bounds[j] then return false; fi;
        od;
        if all_nonneg and all_nonpos then return false; fi;
        if not all_nonneg and not all_nonpos then return false; fi;
        return true;
    end;
    
    # F4 扩展根: a0 = -theta, a1, a2, a3, a4
    F4_ExtendedRoots := [
        [-2, -3, -4, -2],  # a0 = -theta
        [1, 0, 0, 0],      # a1
        [0, 1, 0, 0],      # a2  
        [0, 0, 1, 0],      # a3
        [0, 0, 0, 1]       # a4
    ];
    all_pos_labels := ["a0","a1","a2","a3","a4"];

    BuildSelectedIndices := function(labels)
        local ids, one_label, pos;
        ids := [];
        for one_label in labels do
            pos := Position(all_pos_labels, one_label);
            if pos <> fail then
                Add(ids, pos);
            fi;
        od;
        if Length(ids) = 0 then
            return [1,2,3,4];
        fi;
        return ids;
    end;
    
    MakeKey := function(roots)
        local s, r, j;
        s := "";
        for r in roots do
            for j in [1..Length(r)] do
                if s <> "" then s := Concatenation(s, ","); fi;
                s := Concatenation(s, String(r[j]));
            od;
            s := Concatenation(s, "|");
        od;
        return s;
    end;

    BuildStates := function()
        local all_states;
        queue := [ ShallowCopy(F4_ExtendedRoots) ];
        seen := [ MakeKey(F4_ExtendedRoots) ];
        all_states := [];
        while Length(queue) > 0 do
            state := queue[1];
            Remove(queue, 1);
            Add(all_states, state);
            for k in [1..4] do
                w_roots := [];
                for i in [1..5] do
                    w_roots[i] := F4_Reflection(state[i], k);
                od;
                key := MakeKey(w_roots);
                if Position(seen, key) = fail then
                    Add(seen, key);
                    Add(queue, w_roots);
                fi;
            od;
        od;
        return all_states;
    end;

    BuildRecord := function(one_state, labels, ids, override_colors)
        local roots_out, roots_vecs, all_roots_out, all_colors_out, one_idx, out_i, local_colors, local_black, local_white, c, black_labels, cls, comp_info;
        roots_out := [];
        roots_vecs := [];
        all_roots_out := [];
        all_colors_out := [];
        for one_idx in [1..Length(one_state)] do
            Add(all_roots_out, F4_FormatRoot(one_state[one_idx]));
            if F4_IsRootBlack(one_state[one_idx]) then
                Add(all_colors_out, 1);
            else
                Add(all_colors_out, 0);
            fi;
        od;
        for one_idx in ids do
            Add(roots_vecs, ShallowCopy(one_state[one_idx]));
            Add(roots_out, F4_FormatRoot(one_state[one_idx]));
        od;
        local_colors := [];
        local_black := 0;
        black_labels := [];
        if override_colors <> fail then
            local_colors := ShallowCopy(override_colors);
            for out_i in [1..Length(local_colors)] do
                c := local_colors[out_i];
                if c = 1 then
                    local_black := local_black + 1;
                    Add(black_labels, labels[out_i]);
                fi;
            od;
        else
            for out_i in [1..Length(ids)] do
                one_idx := ids[out_i];
                if F4_IsRootBlack(one_state[one_idx]) then
                    Add(local_colors, 1);
                    local_black := local_black + 1;
                    Add(black_labels, labels[out_i]);
                else
                    Add(local_colors, 0);
                fi;
            od;
        fi;
        local_white := Length(ids) - local_black;
        cls := ValueGlobal("ClassifyRootSystem")(roots_vecs, F4_InnerProduct);
        comp_info := [];
        if PositionSublist(cls.type_str, "+") <> fail then
            comp_info := cls.components;
        fi;
        return rec(
            roots_list := roots_out,
            all_roots_list := all_roots_out,
            all_colors_list := all_colors_out,
            roots_vecs := roots_vecs,
            colors_list := local_colors,
            black_nodes := local_black,
            white_nodes := local_white,
            black_labels := black_labels,
            type := cls.type_str,
            components := comp_info,
            pos_labels := labels,
            all_pos_labels := ShallowCopy(all_pos_labels)
        );
    end;

    results := [];

    if IsRecord(custom_roots) and IsBound(custom_roots.mode) and custom_roots.mode = "SCAN" then
        selected_pos_labels := ["a0","a1","a2","a3"];
        if IsBound(custom_roots.selected_pos_labels) and IsList(custom_roots.selected_pos_labels) and Length(custom_roots.selected_pos_labels) > 0 then
            selected_pos_labels := ShallowCopy(custom_roots.selected_pos_labels);
        fi;
        selected_indices := BuildSelectedIndices(selected_pos_labels);
        target_black_labels := fail;
        if IsBound(custom_roots.target_black_labels) and IsList(custom_roots.target_black_labels) then
            target_black_labels := ShallowCopy(custom_roots.target_black_labels);
        fi;
        required_black_count := fail;
        if IsBound(custom_roots.target_black) and IsInt(custom_roots.target_black) then
            required_black_count := custom_roots.target_black;
        elif IsBound(custom_roots.required_black_count) and IsInt(custom_roots.required_black_count) then
            required_black_count := custom_roots.required_black_count;
        fi;
        target_white_count := fail;
        if IsBound(custom_roots.target_white) and IsInt(custom_roots.target_white) then
            target_white_count := custom_roots.target_white;
        fi;
        target_colors := fail;
        if custom_colors <> fail and IsList(custom_colors) then
            target_colors := ShallowCopy(custom_colors);
        elif IsBound(custom_roots.target_colors) and IsList(custom_roots.target_colors) then
            target_colors := ShallowCopy(custom_roots.target_colors);
        fi;
        target_color_patterns := fail;
        if IsBound(custom_roots.target_color_patterns) and IsList(custom_roots.target_color_patterns) and Length(custom_roots.target_color_patterns) > 0 then
            target_color_patterns := ShallowCopy(custom_roots.target_color_patterns);
        fi;
        states := BuildStates();
        idx := 0;
        for state in states do
            type_info := BuildRecord(state, selected_pos_labels, selected_indices, fail);
            colors_list := type_info.colors_list;
            black_count := type_info.black_nodes;
            white_count := type_info.white_nodes;
            black_labels_actual := type_info.black_labels;
            is_match := true;
            if target_black_labels <> fail then
                if Length(black_labels_actual) <> Length(target_black_labels) then
                    is_match := false;
                else
                    for label in target_black_labels do
                        if Position(black_labels_actual, label) = fail then
                            is_match := false;
                            break;
                        fi;
                    od;
                fi;
            fi;
            if is_match and required_black_count <> fail and black_count <> required_black_count then
                is_match := false;
            fi;
            if is_match and target_white_count <> fail and white_count <> target_white_count then
                is_match := false;
            fi;
            if is_match and target_colors <> fail then
                if Length(target_colors) <> Length(colors_list) or target_colors <> colors_list then
                    is_match := false;
                fi;
            fi;
            if is_match and target_color_patterns <> fail then
                is_match := false;
                for p in target_color_patterns do
                    if IsList(p) and Length(p) = Length(colors_list) and p = colors_list then
                        is_match := true;
                        break;
                    fi;
                od;
            fi;
            if is_match then
                idx := idx + 1;
                Add(results, rec(
                    idx := idx,
                    roots_list := type_info.roots_list,
                    all_roots_list := type_info.all_roots_list,
                    all_colors_list := type_info.all_colors_list,
                    roots_vecs := type_info.roots_vecs,
                    colors_list := colors_list,
                    components := type_info.components,
                    black_nodes := black_count,
                    white_nodes := white_count,
                    black_labels := type_info.black_labels,
                    type := type_info.type,
                    pos_labels := type_info.pos_labels,
                    all_pos_labels := type_info.all_pos_labels
                ));
            fi;
        od;
        return results;
    fi;

    if custom_roots <> fail and IsList(custom_roots) then
        selected_pos_labels := ShallowCopy(custom_roots);
        selected_indices := BuildSelectedIndices(selected_pos_labels);
        if Length(selected_indices) <> Length(selected_pos_labels) then
            Error("F4 自定义模式中存在非法位置标签");
        fi;
        if custom_colors <> fail then
            if not IsList(custom_colors) or Length(custom_colors) <> Length(selected_pos_labels) then
                Error("F4 自定义模式中颜色数量与位置数量不匹配");
            fi;
            target_colors := ShallowCopy(custom_colors);
        else
            target_colors := fail;
        fi;
        state := ShallowCopy(F4_ExtendedRoots);
        type_info := BuildRecord(state, selected_pos_labels, selected_indices, target_colors);
        Add(results, rec(
            idx := 0,
            roots_list := type_info.roots_list,
            all_roots_list := type_info.all_roots_list,
            all_colors_list := type_info.all_colors_list,
            roots_vecs := type_info.roots_vecs,
            colors_list := type_info.colors_list,
            components := type_info.components,
            black_nodes := type_info.black_nodes,
            white_nodes := type_info.white_nodes,
            black_labels := type_info.black_labels,
            type := type_info.type,
            pos_labels := type_info.pos_labels,
            all_pos_labels := type_info.all_pos_labels
        ));
        return results;
    fi;

    return FindF4_B4_Subsystems(
        rec(
            mode := "SCAN",
            selected_pos_labels := ["a0","a1","a2","a3"],
            target_black := 1,
            target_white := 3
        ),
        fail
    );
end;;

# =============================================================================
# B4 根系统生成函数
# =============================================================================

# GenerateB4RootSystem: 生成 B4 的完整根系统
# B4 有 4 个简单根，完整根系统包含 16 个根（8 个正根 + 8 个负根）
# 输入: 4 个简单根（字符串格式）
# 输出: B4 的所有根列表
GenerateB4RootSystem := function(simple_roots)
    local all_roots, i, j, k, l, root, coeff, root_str,
          e1, e2, e3, e4, # 标准基向量
          simple_vecs,    # 简单根的向量表示
          positive_roots, # 正根
          cartan_matrix;
    
    # 如果输入是字符串列表，先转换为向量表示
    simple_vecs := [];
    for i in [1..Length(simple_roots)] do
        # 这里假设 simple_roots 已经是向量形式
        # 如果需要解析字符串，可以添加解析逻辑
        if IsList(simple_roots[i]) and Length(simple_roots[i]) = 4 then
            Add(simple_vecs, simple_roots[i]);
        else
            # 如果是字符串，使用默认的 B4 简单根
            # B4 的标准简单根: e1-e2, e2-e3, e3-e4, e4
            simple_vecs := [
                [1, -1, 0, 0],   # alpha1 = e1-e2
                [0, 1, -1, 0],   # alpha2 = e2-e3  
                [0, 0, 1, -1],   # alpha3 = e3-e4
                [0, 0, 0, 1]     # alpha4 = e4
            ];
            break;
        fi;
    od;
    
    # 生成所有正根
    positive_roots := [];
    
    # 简单根
    for i in [1..4] do
        Add(positive_roots, ShallowCopy(simple_vecs[i]));
    od;
    
    # 高度为 2 的根
    # alpha1 + alpha2
    root := [];
    for j in [1..4] do root[j] := simple_vecs[1][j] + simple_vecs[2][j]; od;
    Add(positive_roots, root);
    
    # alpha2 + alpha3
    root := [];
    for j in [1..4] do root[j] := simple_vecs[2][j] + simple_vecs[3][j]; od;
    Add(positive_roots, root);
    
    # alpha3 + alpha4
    root := [];
    for j in [1..4] do root[j] := simple_vecs[3][j] + simple_vecs[4][j]; od;
    Add(positive_roots, root);
    
    # 高度为 3 的根
    # alpha1 + alpha2 + alpha3
    root := [];
    for j in [1..4] do root[j] := simple_vecs[1][j] + simple_vecs[2][j] + simple_vecs[3][j]; od;
    Add(positive_roots, root);
    
    # alpha2 + alpha3 + alpha4
    root := [];
    for j in [1..4] do root[j] := simple_vecs[2][j] + simple_vecs[3][j] + simple_vecs[4][j]; od;
    Add(positive_roots, root);
    
    # 高度为 4 的根
    # alpha1 + alpha2 + alpha3 + alpha4
    root := [];
    for j in [1..4] do root[j] := simple_vecs[1][j] + simple_vecs[2][j] + simple_vecs[3][j] + simple_vecs[4][j]; od;
    Add(positive_roots, root);
    
    # alpha2 + 2*alpha3 + alpha4 (特殊根)
    root := [];
    for j in [1..4] do root[j] := simple_vecs[2][j] + 2*simple_vecs[3][j] + simple_vecs[4][j]; od;
    Add(positive_roots, root);
    
    # 生成所有根（正根 + 负根）
    all_roots := [];
    
    # 添加正根
    for i in [1..Length(positive_roots)] do
        root_str := "";
        for j in [1..4] do
            coeff := positive_roots[i][j];
            if coeff <> 0 then
                if coeff > 0 and root_str <> "" then root_str := Concatenation(root_str, "+"); fi;
                if coeff = 1 then root_str := Concatenation(root_str, "e", String(j));
                elif coeff = -1 then root_str := Concatenation(root_str, "-e", String(j));
                else root_str := Concatenation(root_str, String(coeff), "e", String(j)); fi;
            fi;
        od;
        if root_str = "" then root_str := "0"; fi;
        Add(all_roots, root_str);
    od;
    
    # 添加负根
    for i in [1..Length(positive_roots)] do
        root_str := "";
        for j in [1..4] do
            coeff := -positive_roots[i][j];
            if coeff <> 0 then
                if coeff > 0 and root_str <> "" then root_str := Concatenation(root_str, "+"); fi;
                if coeff = 1 then root_str := Concatenation(root_str, "e", String(j));
                elif coeff = -1 then root_str := Concatenation(root_str, "-e", String(j));
                else root_str := Concatenation(root_str, String(coeff), "e", String(j)); fi;
            fi;
        od;
        if root_str = "" then root_str := "0"; fi;
        Add(all_roots, root_str);
    od;
    
    return all_roots;
end;;
# 简单反射
# s_alpha(v) = v - 2 (v,alpha)/(alpha,alpha) * alpha
# 对 G2 根系该系数为整数，因此反射后仍落在根格中
Reflection := function(v, alpha)
    local coeff;
    coeff := 2 * InnerProduct(v, alpha) / InnerProduct(alpha, alpha);
    return v - coeff * alpha;
end;;

# 2. Vogan 图着色规则 (G2 Split Real Form)
# 规则：以第二坐标的奇偶性作为黑/白判定（黑=1，白=0）
IsRootBlack := function(root)
    local n1;
    n1 := root[1];
    if n1 < 0 then n1 := -n1; fi;
    return (n1 mod 2) <> 0;
end;;

# 3. 格式化输出辅助函数
# 将 [x,y] 形式的根向量格式化为 "xa1+ya2" 字符串（省略 1，保留 ±）
FormatRoot := function(v)
    local s;
    s := "";
    if v[1] <> 0 then 
        if v[1] = 1 then s := Concatenation(s, "a1"); 
        elif v[1] = -1 then s := Concatenation(s, "-a1");
        else s := Concatenation(s, String(v[1]), "a1"); fi;
    fi;
    if v[2] <> 0 then
        if v[1] <> 0 and v[2] > 0 then s := Concatenation(s, "+"); fi;
        if v[2] = 1 then s := Concatenation(s, "a2");
        elif v[2] = -1 then s := Concatenation(s, "-a2");
        else s := Concatenation(s, String(v[2]), "a2"); fi;
    fi;
    if s = "" then s := "0"; fi;
    return s;
end;;

# =============================================================================
# 第二部分: 核心分析函数
# =============================================================================

# 返回符合条件的子系统列表
# 输入:
#   custom_roots:
#     - 记录且 mode="SCAN": 使用扫描模式 (根据模板与目标黑白数筛选)
#     - 列表: 指定根字符串列表（如 ["a1", "-3a1-2a2"]）
#     - fail: 使用默认扩展根 (alpha2, -theta) 扫描
#   custom_colors:
#     - 与 custom_roots 为列表时可给出颜色列表（1=黑,0=白），或传 fail 自动推断
# 输出:
#   子系统记录列表，每个记录包含:
#     idx, roots_list, colors_list, components, black_nodes, white_nodes, type
FindG2Subsystems := function(custom_roots, custom_colors)
    local WeylWords, idx, word, w_alpha1, w_alpha2, k, ref_alpha, w_theta, neg_w_theta,
          is_alpha2_black, is_ext_black, subsystem, type_str, results,
          black_nodes, white_nodes, root_strs, root_vecs, i, color, root_v,
          type_info, colors_list, c1, c2,
          target_black, target_white,
          all_roots, idx_count, r1, r2, ip, j,
          desired_type, desired_long_long, t1, t2, len1, len2, seen, key,
          l1, l2, tmp, tmpc;

    # 分类器和 K-Orbit 分类器应该在主流程中已经加载
    
    results := [];
    
    if IsRecord(custom_roots) and IsBound(custom_roots.mode) and custom_roots.mode = "SCAN" then
        # 扫描模式: 根据实验意图扫描子系统
        target_black := custom_roots.target_black;
        target_white := custom_roots.target_white;
        
        # 推断期望的子系统类型
        desired_type := fail;           # "A2" | "A1 + A1" | fail
        desired_long_long := false;     # 是否要求两根均为长根
        if IsRecord(custom_roots) and IsBound(custom_roots.template) and Length(custom_roots.template) = 2 then
            t1 := custom_roots.template[1];
            t2 := custom_roots.template[2];
            len1 := InnerProduct(t1, t1);
            len2 := InnerProduct(t2, t2);
            if target_black = 1 and target_white = 1 and len1 = 6 and len2 = 6 then
                # 长根 + 长根, 一黑一白 -> A2
                desired_type := "A2";
                desired_long_long := true;
            elif target_black = 2 and target_white = 0 then
                # 两黑 -> A1 + A1（正交）
                desired_type := "A1 + A1";
            fi;
        fi;
        
        if desired_type = "A2" and desired_long_long then
            WeylWords := [
                [], [1], [2], [1,2], [2,1], [1,2,1], [2,1,2], [1,2,1,2], [2,1,2,1], [1,2,1,2,1], [2,1,2,1,2], [1,2,1,2,1,2]
            ];
            results := [];
            idx_count := 0;
            seen := [];
            for idx in [1..Length(WeylWords)] do
                word := WeylWords[idx];
                w_alpha1 := [1, 0];
                w_alpha2 := [0, 1];
                for k in Reversed(word) do
                    ref_alpha := SimpleRoots[k];
                    w_alpha1 := Reflection(w_alpha1, ref_alpha);
                    w_alpha2 := Reflection(w_alpha2, ref_alpha);
                od;
                w_theta := 3 * w_alpha1 + 2 * w_alpha2;
                neg_w_theta := -w_theta;
                
                c1 := 0; if IsRootBlack(w_alpha2) then c1 := 1; fi;
                c2 := 0; if IsRootBlack(neg_w_theta) then c2 := 1; fi;
                if (c1 + c2) <> target_black or ((1-c1) + (1-c2)) <> target_white then
                    continue;
                fi;
                
                type_info := ValueGlobal("ClassifyRootSystem")([w_alpha2, neg_w_theta], InnerProduct);
                if type_info.type_str <> desired_type then
                    continue;
                fi;
                
                r1 := neg_w_theta;
                r2 := w_alpha2;
                root_strs := [ FormatRoot(r1), FormatRoot(r2) ];
                key := Concatenation(root_strs[1], "|", root_strs[2]);
                if Position(seen, key) <> fail then
                    continue;
                fi;
                Add(seen, key);
                idx_count := idx_count + 1;
                colors_list := [ c2, c1 ];
                
                Add(results, rec(
                    idx := idx_count,
                    roots_list := root_strs,
                    colors_list := colors_list,
                    components := type_info.components,
                    black_nodes := (c1 + c2),
                    white_nodes := ((1-c1) + (1-c2)),
                    type := type_info.type_str,
                    w_alpha2 := FormatRoot(w_alpha2),
                    neg_w_theta := FormatRoot(neg_w_theta)
                ));
            od;
            return results;
        fi;
        
        # 1. 生成 G2 所有 12 个根
        all_roots := [
            [1, 0], [0, 1], [1, 1], [2, 1], [3, 1], [3, 2],
            [-1, 0], [0, -1], [-1, -1], [-2, -1], [-3, -1], [-3, -2]
        ];
        
        results := [];
        idx_count := 0;
        
        # 2. 遍历所有唯一的根对 (r1, r2)
        for i in [1..Length(all_roots)] do
            for j in [i+1..Length(all_roots)] do
                r1 := all_roots[i];
                r2 := all_roots[j];
                
                # 计算颜色 (1=黑,0=白)
                c1 := 0; if IsRootBlack(r1) then c1 := 1; fi;
                c2 := 0; if IsRootBlack(r2) then c2 := 1; fi;
                
                if (c1 + c2) <> target_black or ((1-c1) + (1-c2)) <> target_white then
                    continue;
                fi;
                
                # 分类系统类型
                type_info := ValueGlobal("ClassifyRootSystem")([r1, r2], InnerProduct);
                
                # 根据期望类型进行过滤
                if desired_type <> fail then
                    if type_info.type_str <> desired_type then
                        continue;
                    fi;
                    # 若期望长-长，进一步过滤长度
                    if desired_long_long then
                        if InnerProduct(r1, r1) <> 6 or InnerProduct(r2, r2) <> 6 then
                            continue;
                        fi;
                    fi;
                fi;
                
                # 对于显式 A1 + A1，可确保正交性
                if desired_type = "A1 + A1" then
                    ip := InnerProduct(r1, r2);
                    if ip <> 0 then
                        continue;
                    fi;
                fi;
                
                l1 := InnerProduct(r1, r1);
                l2 := InnerProduct(r2, r2);
                if l1 < l2 then
                    tmp := r1; r1 := r2; r2 := tmp;
                    tmpc := c1; c1 := c2; c2 := tmpc;
                fi;
                idx_count := idx_count + 1;
                root_strs := [ FormatRoot(r1), FormatRoot(r2) ];
                colors_list := [ c1, c2 ];
                
                Add(results, rec(
                    idx := idx_count,
                    roots_list := root_strs,
                    colors_list := colors_list,
                    components := type_info.components,
                    black_nodes := (c1 + c2),
                    white_nodes := ((1-c1) + (1-c2)),
                    type := type_info.type_str,
                    # 兼容字段
                    w_alpha2 := root_strs[1],
                    neg_w_theta := root_strs[2]
                ));
            od;
        od;
        return results;
    fi;
    
    # 模式 1: 自定义根集模式
    if custom_roots <> fail then
        # 解析根向量
        root_vecs := [];
        for i in [1..Length(custom_roots)] do
            Add(root_vecs, ValueGlobal("ParseRootString")(custom_roots[i]));
        od;
        
        # 确定颜色
        black_nodes := 0;
        white_nodes := 0;
        
        if custom_colors <> fail then
            # 使用用户指定的颜色
            if Length(custom_colors) <> Length(custom_roots) then
                Error("颜色数量与根数量不匹配");
            fi;
            for color in custom_colors do
                if color = 1 then black_nodes := black_nodes + 1;
                else white_nodes := white_nodes + 1;
                fi;
            od;
        else
            # 使用默认 G2 Split 规则推断颜色
            for root_v in root_vecs do
                if IsRootBlack(root_v) then black_nodes := black_nodes + 1;
                else white_nodes := white_nodes + 1;
                fi;
            od;
        fi;
        
        # 判定类型
        type_info := ValueGlobal("ClassifyRootSystem")(root_vecs, InnerProduct);
        type_str := type_info.type_str;
        
        # 构造 colors_list 如果用户未提供
        if custom_colors = fail then
            colors_list := [];
            for root_v in root_vecs do
                if IsRootBlack(root_v) then 
                    Add(colors_list, 1);
                else 
                    Add(colors_list, 0);
                fi;
            od;
        else
            colors_list := custom_colors;
        fi;
        
        # 构造单个结果
        Add(results, rec(
            idx := 0, # 自定义模式无 Weyl 索引
            w_alpha2 := custom_roots[1], # 假设至少有一个根，用于显示
            neg_w_theta := custom_roots[Length(custom_roots)], # 显示最后一个
            
            # 扩展字段: roots_list
            roots_list := custom_roots,
            
            colors_list := colors_list,
            components := type_info.components,
            
            black_nodes := black_nodes,
            white_nodes := white_nodes,
            type := type_str
        ));
        
        return results;
    fi;

    # 模式 2: 默认 G2 扩展根搜索 (原逻辑)
    # 对 (alpha2, -theta) 的全 Weyl 轨道进行扫描，筛选一黑一白的配置
    # 生成 Weyl 群 (12 个元素)
    WeylWords := [
        [], [1], [2], [1,2], [2,1], [1,2,1], [2,1,2], [1,2,1,2], [2,1,2,1], [1,2,1,2,1], [2,1,2,1,2], [1,2,1,2,1,2]
    ];
    
    for idx in [1..Length(WeylWords)] do
        word := WeylWords[idx];
        
        # 初始基
        w_alpha1 := [1, 0];
        w_alpha2 := [0, 1];
        
        # 应用 Weyl 变换 (从右向左)
        for k in Reversed(word) do
            ref_alpha := SimpleRoots[k];
            w_alpha1 := Reflection(w_alpha1, ref_alpha);
            w_alpha2 := Reflection(w_alpha2, ref_alpha);
        od;
        
        # 计算 w(theta) 和 -w(theta)
        w_theta := 3 * w_alpha1 + 2 * w_alpha2;
        neg_w_theta := -w_theta;
        
        # 颜色判断
        is_alpha2_black := IsRootBlack(w_alpha2);
        is_ext_black := IsRootBlack(neg_w_theta);
        
        # 筛选: 一黑一白
        if is_alpha2_black <> is_ext_black then
            # 统计颜色
            black_nodes := 0;
            white_nodes := 0;
            if is_alpha2_black then black_nodes := black_nodes + 1; else white_nodes := white_nodes + 1; fi;
            if is_ext_black then black_nodes := black_nodes + 1; else white_nodes := white_nodes + 1; fi;
            
            # 子系统类型判断
            subsystem := [ w_alpha2, neg_w_theta ];
            type_info := ValueGlobal("ClassifyRootSystem")(subsystem, InnerProduct);
            
            c1 := 0; if is_alpha2_black then c1 := 1; fi;
            c2 := 0; if is_ext_black then c2 := 1; fi;
            
            Add(results, rec(
                idx := idx,
                w_alpha2 := FormatRoot(w_alpha2),
                neg_w_theta := FormatRoot(neg_w_theta),
                roots_list := [ FormatRoot(w_alpha2), FormatRoot(neg_w_theta) ],
                colors_list := [ c1, c2 ],
                components := type_info.components,
                black_nodes := black_nodes,
                white_nodes := white_nodes,
                type := type_info.type_str
            ));
        fi;
    od;
    
    return results;
end;;
