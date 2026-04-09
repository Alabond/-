#############################################################################
# 文件: Main_Pipeline.g
#############################################################################

Print("BOOT\n");

# 允许外部覆盖输出文件名
if not IsBound(GlobalOutputDir) then
    GlobalOutputDir := "outputs";
fi;
out_file := Concatenation(GlobalOutputDir, "/F4_pipeline_output.txt");
if IsBound(GlobalOutFile) then
    out_file := GlobalOutFile;
fi;
PrintTo(out_file, "");
LogTo(out_file);

if not IsBound(GlobalAmbientType) then
    GlobalAmbientType := "F4";
fi;
if not IsBound(GlobalEnableLogFile) then
    GlobalEnableLogFile := false;
fi;
if not IsBound(GlobalF4CompositeOutputStyle) and GlobalAmbientType = "F4" then
    GlobalF4CompositeOutputStyle := "single_like";
fi;

# 首先加载基础工具模块
Print("LOAD 02_Subsystem_Classification\n");
Read("../02_Subsystem_Classification/Subsystem_Classifier.g");
Print("LOAD 04_Orbit_Analysis/K_Orbit_Classifier\n");
Read("../04_Orbit_Analysis/K_Orbit_Classifier.g");

# 然后加载其他模块
Print("LOAD 01_Subsystem_Generation\n");
if GlobalEnableLogFile and not IsBound(GlobalLogfileName) then
    GlobalLogfileName := Concatenation(GlobalOutputDir, "/pipeline.log");
fi;
Read("../01_Subsystem_Generation/G2_Extended_Subsystem_Lib.g");
Print("LOAD 04_Orbit_Analysis/Noticed_Orbits\n");
Read("../04_Orbit_Analysis/Noticed_Orbits.g");
Print("LOAD 03_Real_Form_Inference\n");
Read("../03_Real_Form_Inference/Real_Form_Inference.g");
Print("LOAD 05_SL2_Triple_Construction\n");
Read("../05_SL2_Triple_Construction/SL2_Triple_Builder.g");

if GlobalEnableLogFile and IsBound(GlobalLogfileName) and GlobalLogfileName <> fail then
    logfile := GlobalLogfileName;
    PrintTo(logfile, "");
else
    logfile := fail;
fi;
AppendPipelineLog := function(text)
    if logfile <> fail then
        AppendTo(logfile, text);
    fi;
end;;
Print("START\n");
AppendPipelineLog("START\n");
EmitSummary := function(key, value)
    Print("SUMMARY ", key, ": ", value, "\n");
    AppendPipelineLog(Concatenation("SUMMARY ", key, " ", String(value), "\n"));
end;;
GetGlobalOrDefault := function(global_name, default_value)
    if IsBoundGlobal(global_name) then
        return ValueGlobal(global_name);
    fi;
    return default_value;
end;;
EmitSummaryIfGlobal := function(key, global_name)
    if IsBoundGlobal(global_name) then
        EmitSummary(key, ValueGlobal(global_name));
    fi;
end;;
F4PipelineInnerProduct := function(v1, v2)
    local M, a, b, res;
    M := [
        [4, -2, 0, 0],
        [-2, 4, -2, 0],
        [0, -2, 2, -1],
        [0, 0, -1, 2]
    ];
    res := 0;
    for a in [1..4] do
        for b in [1..4] do
            res := res + v1[a] * M[a][b] * v2[b];
        od;
    od;
    return res;
end;;
ResolveBCTypeByRootLengths := function(component_vecs)
    local lens, min_len, max_len, short_count, long_count;
    if not IsList(component_vecs) or Length(component_vecs) < 2 then
        return fail;
    fi;
    lens := List(component_vecs, v -> F4PipelineInnerProduct(v, v));
    min_len := Minimum(lens);
    max_len := Maximum(lens);
    if min_len = max_len then
        return fail;
    fi;
    short_count := Length(Filtered(lens, x -> x = min_len));
    long_count := Length(Filtered(lens, x -> x = max_len));
    if long_count = 1 and short_count = Length(lens) - 1 then
        return "C";
    fi;
    if short_count = 1 and long_count = Length(lens) - 1 then
        return "B";
    fi;
    return fail;
end;;

ResolveComponentTypeByRootData := function(raw_type, component_vecs)
    local slash_pos, resolved_bc_type, rank_part;
    slash_pos := Position(raw_type, '/');
    if slash_pos <> fail and slash_pos > 1 then
        if IsList(component_vecs) and Length(component_vecs) > 0 then
            resolved_bc_type := ResolveBCTypeByRootLengths(component_vecs);
            if resolved_bc_type <> fail then
                rank_part := raw_type{[2..slash_pos-1]};
                return Concatenation(resolved_bc_type, rank_part);
            fi;
        fi;
        return raw_type{[1..slash_pos-1]};
    fi;
    return raw_type;
end;;

ResolveSubsystemDisplayType := function(s)
    local pieces, comp, comp_vecs, idx, resolved_type, result, i;
    if not (IsBound(s.components) and IsList(s.components) and Length(s.components) > 0 and IsBound(s.roots_vecs)) then
        return s.type;
    fi;
    pieces := [];
    for comp in s.components do
        comp_vecs := [];
        for idx in comp.indices do
            Add(comp_vecs, s.roots_vecs[idx]);
        od;
        resolved_type := ResolveComponentTypeByRootData(comp.type, comp_vecs);
        Add(pieces, resolved_type);
    od;
    Sort(pieces);
    result := "";
    for i in [1..Length(pieces)] do
        if i > 1 then
            result := Concatenation(result, " + ");
        fi;
        result := Concatenation(result, pieces[i]);
    od;
    return result;
end;;

InferRealForm_A3_ByPosition := function(component_vecs, component_colors)
    local C, i, j, deg, middle_idx, blacks, whites, minority_idx;
    if Length(component_vecs) <> 3 or Length(component_colors) <> 3 then
        return fail;
    fi;
    C := UC_ComputeCartanMatrix(component_vecs, F4PipelineInnerProduct);
    deg := [];
    for i in [1..3] do
        Add(deg, 0);
        for j in [1..3] do
            if i <> j and (C[i][j] <> 0 or C[j][i] <> 0) then
                deg[i] := deg[i] + 1;
            fi;
        od;
    od;
    middle_idx := fail;
    for i in [1..3] do
        if deg[i] = 2 then
            middle_idx := i;
            break;
        fi;
    od;
    if middle_idx = fail then
        return fail;
    fi;
    blacks := Filtered([1..3], i -> component_colors[i] = 1);
    whites := Filtered([1..3], i -> component_colors[i] = 0);
    if Length(blacks) = 0 or Length(whites) = 0 then
        return fail;
    fi;
    if Length(blacks) <= Length(whites) then
        minority_idx := blacks[1];
    else
        minority_idx := whites[1];
    fi;
    if minority_idx = middle_idx then
        return rec(type := "su_pq", params := [2, 2], desc := "su(2,2)");
    fi;
    return rec(type := "su_pq", params := [3, 1], desc := "su(3,1)");
end;;

InferRealForm_C3_ByPosition := function(component_vecs, component_colors)
    local lens, max_len, long_idxs, black_idxs;
    if Length(component_vecs) <> 3 or Length(component_colors) <> 3 then
        return fail;
    fi;
    lens := List(component_vecs, v -> F4PipelineInnerProduct(v, v));
    max_len := Maximum(lens);
    long_idxs := Filtered([1..3], i -> lens[i] = max_len);
    if Length(long_idxs) <> 1 then
        return fail;
    fi;
    black_idxs := Filtered([1..3], i -> component_colors[i] = 1);
    if Length(black_idxs) = 1 and black_idxs[1] = long_idxs[1] then
        return rec(type := "sp_R", params := [3], desc := "sp(3, R)");
    fi;
    return fail;
end;;

InferRealForm_C2_ByPosition := function(component_vecs, component_colors)
    local lens, max_len, long_idxs, black_idxs;
    if Length(component_vecs) <> 2 or Length(component_colors) <> 2 then
        return fail;
    fi;
    lens := List(component_vecs, v -> F4PipelineInnerProduct(v, v));
    max_len := Maximum(lens);
    long_idxs := Filtered([1..2], i -> lens[i] = max_len);
    if Length(long_idxs) <> 1 then
        return fail;
    fi;
    black_idxs := Filtered([1..2], i -> component_colors[i] = 1);
    if Length(black_idxs) = 1 and black_idxs[1] = long_idxs[1] then
        return rec(type := "sp_R", params := [2], desc := "sp(2, R)");
    fi;
    return fail;
end;;
GlobalNormalTripleTotalChecks := 0;
GlobalNormalTriplePassChecks := 0;
GlobalF4LabelStats := [];
GlobalF4LabelFallbackCount := 0;
GlobalF4LabelFallbackMaxErr := 3;

# 2. 步骤 1: G2 子系统筛选与配置
Print("\n==================================================\n");
Print("步骤 1: 子系统选取与配置\n");
Print("==================================================\n");
# 允许外部覆盖简介文本
if IsBound(GlobalIntroText) then
    Print(GlobalIntroText, "\n");
else
    Print("研究对象: 例外实李代数 F4(4)\n");
    Print("扩展示验: 左边含 -theta 的 B4\n");
fi;

# 用户配置区域（允许外部注入实验；默认仅运行 F4_B4 实验）
if IsBound(GlobalExperiments) then
    experiments := GlobalExperiments;
else
    experiments := [
        rec(
            desc := "Experiment: F4(4) 选取[a0,a1,a2,a4] 且 a1,a4 涂黑",
            roots := rec(
                mode := "SCAN",
                selected_pos_labels := ["a0","a1","a2","a4"],
                target_black_labels := ["a1","a4"]
            ),
            colors := fail
        )
    ];
fi;

all_subsystems := [];

for exp in experiments do
    Print("\n>>> 运行实验: ", exp.desc, "\n");
    if IsRecord(exp.roots) and exp.roots.mode = "SCAN" then
        Print("  Mode: Scan\n");
        if IsBound(exp.roots.template) then
            Print("  Template: ", exp.roots.template, "\n");
        fi;
        if IsBound(exp.roots.target_black) and IsBound(exp.roots.target_white) then
            Print("  Target: ", exp.roots.target_black, "B ", exp.roots.target_white, "W\n");
        fi;
        if IsBound(exp.roots.selected_pos_labels) then
            Print("  Selected: ", exp.roots.selected_pos_labels, "\n");
        fi;
        if IsBound(exp.roots.target_black_labels) then
            Print("  Target black labels: ", exp.roots.target_black_labels, "\n");
        fi;
        if IsBound(exp.roots.target_colors) then
            Print("  Target colors: ", exp.roots.target_colors, "\n");
        elif IsBound(exp.roots.target_color_patterns) then
            Print("  Target color patterns: ", exp.roots.target_color_patterns, "\n");
        elif IsBound(exp.colors) and exp.colors <> fail then
            Print("  Target colors: ", exp.colors, "\n");
        fi;
    elif exp.roots <> fail then
        Print("  Roots: ", exp.roots, "\n");
        Print("  Colors: ", exp.colors, "\n");
    else
        Print("  Mode: Auto Scan\n");
    fi;
    
    if GlobalAmbientType = "F4" then
        if IsBound(exp.roots) then
            subs := FindF4_B4_Subsystems(exp.roots, exp.colors);
        else
            subs := FindF4_B4_Subsystems(fail, fail);
        fi;
    else
        subs := FindG2Subsystems(exp.roots, exp.colors);
    fi;
    Print("  发现 ", Length(subs), " 个子系统配置。\n");
    AppendPipelineLog(Concatenation("EXP: ", exp.desc, ", subs=", String(Length(subs)), "\n"));
    
    # 为子系统添加实验描述标记，以便区分
    for s in subs do
        s.exp_desc := exp.desc;
    od;
    
    Append(all_subsystems, subs);
od;

subsystems := all_subsystems;

# 文本对齐辅助（右侧填充空格至固定宽度）
PadRight := function(s, w)
    local t, len_s;
    t := s;
    len_s := Length(s);
    
    if len_s >= w then return s; fi;
    
    while Length(t) < w do t := Concatenation(t, " "); od;
    return t;
end;;

# 打印筛选结果概览
Print("\n筛选结果概览:\n");
header := Concatenation(
    PadRight("索引", 6), "| ",
    PadRight("Roots", 34), "| ", 
    PadRight("颜色(B/W)", 20), "| ",  # 增加宽度以容纳位置信息
    PadRight("类型", 10)
);
Print(header, "\n");
Print("--------------------------------------------------------------------------\n");

for s in subsystems do
    color_desc := Concatenation(String(s.black_nodes), "B / ", String(s.white_nodes), "W");
    
    # 格式化根列表显示
    roots_display := "";
    if IsBound(s.roots_list) then
        roots_display := String(s.roots_list);
    else
        roots_display := Concatenation(s.w_alpha2, ", ", s.neg_w_theta);
    fi;
    
    # 添加涂黑位置详细信息
        black_pos_desc := "";
        if IsBound(s.colors_list) then
            black_positions := [];
            root_labels := [];
            if IsBound(s.pos_labels) then
                root_labels := s.pos_labels;
            elif IsBound(s.roots_list) then
                root_labels := s.roots_list;
            fi;
            for i in [1..Length(s.colors_list)] do
                if s.colors_list[i] = 1 then
                    if i <= Length(root_labels) then
                        Add(black_positions, root_labels[i]);
                    fi;
                fi;
            od;
            if Length(black_positions) > 0 then
                pos_str := "";
                for j in [1..Length(black_positions)] do
                    if j > 1 then
                        pos_str := Concatenation(pos_str, ",");
                    fi;
                    pos_str := Concatenation(pos_str, black_positions[j]);
                od;
                black_pos_desc := Concatenation(" [", pos_str, "]");
            fi;
        fi;
    
    display_type := ResolveSubsystemDisplayType(s);
    Print(
        PadRight(String(s.idx), 6), "| ",
        PadRight(roots_display, 34), "| ",
        PadRight(Concatenation(color_desc, black_pos_desc), 20), "| ",
        PadRight(display_type, 10), "\n"
    );
od;

# 3. 步骤 2: 自动分析每个唯一子系统类型的 Noticed Orbits
Print("\n==================================================\n");
Print("步骤 2: 自动分析 Noticed Orbits\n");
Print("==================================================\n");

# 辅助函数: 计算列表族的笛卡尔积
# 例如 [[a,b],[1,2]] -> [[a,1],[a,2],[b,1],[b,2]]
CartesianProductList := function(list_of_lists)
    local result, list, item, new_result, temp;
    if Length(list_of_lists) = 0 then return [[]]; fi;
    
    result := [[]];
    for list in list_of_lists do
        new_result := [];
        for temp in result do
            for item in list do
                Add(new_result, Concatenation(temp, [item]));
            od;
        od;
        result := new_result;
    od;
    return result;
end;;

F4PipelineFormatRoot := function(v)
    local s, i, coeff, abs_coeff, first;
    s := "";
    first := true;
    for i in [1..4] do
        coeff := v[i];
        if coeff <> 0 then
            if not first and coeff > 0 then
                s := Concatenation(s, "+");
            fi;
            abs_coeff := coeff;
            if abs_coeff < 0 then
                abs_coeff := -abs_coeff;
            fi;
            if coeff = -1 then
                s := Concatenation(s, "-a", String(i));
            elif coeff = 1 then
                s := Concatenation(s, "a", String(i));
            elif coeff < 0 then
                s := Concatenation(s, "-", String(abs_coeff), "a", String(i));
            else
                s := Concatenation(s, String(coeff), "a", String(i));
            fi;
            first := false;
        fi;
    od;
    if s = "" then
        s := "0";
    fi;
    return s;
end;;

F4ReflectVectorBySimpleIndex := function(v, idx)
    local beta, beta_len2, coeff, i, out;
    if idx = 0 then
        beta := [-2, -3, -4, -2];
    else
        beta := [0, 0, 0, 0];
        beta[idx] := 1;
    fi;
    beta_len2 := F4PipelineInnerProduct(beta, beta);
    coeff := 2 * F4PipelineInnerProduct(v, beta) / beta_len2;
    out := [];
    for i in [1..4] do
        Add(out, v[i] - coeff * beta[i]);
    od;
    return out;
end;;

F4ReflectVectorByRoot := function(v, beta)
    local beta_len2, coeff, i, out;
    beta_len2 := F4PipelineInnerProduct(beta, beta);
    coeff := 2 * F4PipelineInnerProduct(v, beta) / beta_len2;
    out := [];
    for i in [1..4] do
        Add(out, v[i] - coeff * beta[i]);
    od;
    return out;
end;;

GetSubsystemRootVectorsF4 := function(s, subsystem_roots)
    if IsBound(s.roots_vecs) and IsList(s.roots_vecs) and Length(s.roots_vecs) > 0 then
        return List(s.roots_vecs, v -> ShallowCopy(v));
    fi;
    if IsBoundGlobal("ParseRootStringFlexible") then
        return List(subsystem_roots, r -> ValueGlobal("ParseRootStringFlexible")(r));
    fi;
    return fail;
end;;

GenerateDeltaJRootSetF4 := function(simple_roots)
    local roots, queue, seen, beta, root, next_root, key_of;
    key_of := v -> F4PipelineFormatRoot(v);
    if simple_roots = fail or not IsList(simple_roots) or Length(simple_roots) = 0 then
        return [];
    fi;
    roots := [];
    seen := [];
    for beta in simple_roots do
        if Position(seen, key_of(beta)) = fail then
            Add(roots, ShallowCopy(beta));
            Add(seen, key_of(beta));
        fi;
    od;
    queue := ShallowCopy(roots);
    while Length(queue) > 0 do
        root := queue[1];
        Remove(queue, 1);
        for beta in simple_roots do
            next_root := F4ReflectVectorByRoot(root, beta);
            if Position(seen, key_of(next_root)) = fail then
                Add(roots, next_root);
                Add(queue, next_root);
                Add(seen, key_of(next_root));
            fi;
        od;
    od;
    return roots;
end;;

CanonicalDeltaJWKOrbitKeyF4 := function(type_str, roots_vecs)
    local generators, queue, seen, best_key, state, state_key, gen, next_state, next_key, delta_j_roots;
    generators := [0, 2, 3, 4];
    state_key := function(vecs)
        local entries, entry, key, sep, v;
        entries := [];
        for v in vecs do
            Add(entries, F4PipelineFormatRoot(v));
        od;
        entries := SortedList(entries);
        key := Concatenation(type_str, "|", String(Length(entries)), "|");
        sep := "";
        for entry in entries do
            key := Concatenation(key, sep, entry);
            sep := ";";
        od;
        return key;
    end;
    if roots_vecs = fail or not IsList(roots_vecs) or Length(roots_vecs) = 0 then
        return Concatenation(type_str, "|0|");
    fi;
    delta_j_roots := GenerateDeltaJRootSetF4(roots_vecs);
    queue := [List(delta_j_roots, v -> ShallowCopy(v))];
    seen := [state_key(queue[1])];
    best_key := seen[1];
    while Length(queue) > 0 do
        state := queue[1];
        Remove(queue, 1);
        for gen in generators do
            next_state := List(state, v -> F4ReflectVectorBySimpleIndex(v, gen));
            next_key := state_key(next_state);
            if Position(seen, next_key) = fail then
                Add(seen, next_key);
                Add(queue, next_state);
                if next_key < best_key then
                    best_key := next_key;
                fi;
            fi;
        od;
    od;
    return best_key;
end;;

GlobalDeltaJWKOrbitSeen := [];
GlobalDeltaJWKOrbitExperimentTag := fail;

RegisterDeltaJWKOrbitRepresentativeF4 := function(subsystem_index, type_str, roots_vecs)
    local key, item;
    if not IsBoundGlobal("GlobalDeltaJWKOrbitSeen") then
        GlobalDeltaJWKOrbitSeen := [];
    fi;
    key := CanonicalDeltaJWKOrbitKeyF4(type_str, roots_vecs);
    for item in GlobalDeltaJWKOrbitSeen do
        if item.key = key then
            return rec(is_new := false, first_index := item.first_index, key := key);
        fi;
    od;
    Add(GlobalDeltaJWKOrbitSeen, rec(key := key, first_index := subsystem_index));
    return rec(is_new := true, first_index := subsystem_index, key := key);
end;;

for s in subsystems do
    if IsBound(s.exp_desc) and GlobalDeltaJWKOrbitExperimentTag <> s.exp_desc then
        GlobalDeltaJWKOrbitSeen := [];
        GlobalDeltaJWKOrbitExperimentTag := s.exp_desc;
    fi;
    Print("\n>>> 模块 1 输出: 子系统分析 (索引 ", s.idx, ")\n");
    if IsBound(s.exp_desc) then
        Print("    实验: ", s.exp_desc, "\n");
    fi;
    display_type := ResolveSubsystemDisplayType(s);
    Print("    类型: ", display_type, "\n");
    
    if IsBound(s.roots_list) then
        Print("    根基: ", s.roots_list, "\n");
        subsystem_roots := s.roots_list;
    else
        Print("    根基: ", s.w_alpha2, ", ", s.neg_w_theta, "\n");
        subsystem_roots := [ s.w_alpha2, s.neg_w_theta ];
    fi;
    
    Print("    颜色配置: ", s.black_nodes, " 黑, ", s.white_nodes, " 白\n");
    
    # 检查是否包含组件信息 (复合子系统)
    if IsBound(s.components) and IsList(s.components) and Length(s.components) > 0 then
        composite_single_like := (IsBoundGlobal("GlobalF4CompositeOutputStyle") and ValueGlobal("GlobalF4CompositeOutputStyle") = "single_like");
        if not composite_single_like then
            Print("    >>> 检测到复合子系统，包含 ", Length(s.components), " 个组件\n");
        fi;
        
        # 1. 对每个组件分别进行实形式推断和 Noticed Orbit 分析
        comp_orbits_list := [];
        comp_roots_list := [];
        comp_infos_list := [];
        
        all_comps_valid := true;
        
        for comp in s.components do
            # comp.indices: 组件在 subsystem_roots 中的索引列表
            # 提取该组件的根
            comp_roots := [];
            comp_vecs := [];
            for idx in comp.indices do
                Add(comp_roots, subsystem_roots[idx]);
                if IsBound(s.roots_vecs) then
                    Add(comp_vecs, s.roots_vecs[idx]);
                fi;
            od;
            Add(comp_roots_list, comp_roots);
            
            # 统计该组件的颜色
            # s.colors_list 包含所有根的颜色 (1=B, 0=W)
            comp_black := 0;
            comp_white := 0;
            if IsBound(s.colors_list) then
                for idx in comp.indices do
                    if s.colors_list[idx] = 1 then comp_black := comp_black + 1;
                    else comp_white := comp_white + 1;
                    fi;
                od;
            else
                # Fallback: 如果没有 colors_list (旧数据结构)，假设是 G2 扩展根 (2个根)
                # 这种情况下 components 应该只有一个? 或者我们需要根据 IsRootBlack 重新计算
                # 这里简单起见，如果 s.colors_list 不存在，我们可能无法处理复合类型
                Print("    警告: 缺少颜色列表，无法处理复合类型组件颜色。\n");
                all_comps_valid := false;
                break;
            fi;
            
            # 解析组件类型（如 "A1","A2"）
            raw_type := ResolveComponentTypeByRootData(comp.type, comp_vecs);
            type_char := raw_type{[1]};
            rank_str := raw_type{[2..Length(raw_type)]};
            rank := Int(rank_str);
            
            if not composite_single_like then
                Print("    --- 组件: ", comp.type, " (根: ", comp_roots, ", 颜色: ", comp_black, "B ", comp_white, "W) ---\n");
            fi;
            
            # 推断实形式
            info := InferRealForm(type_char, rank, comp_black, comp_white);
            comp_colors := [];
            if IsBound(s.colors_list) then
                for idx in comp.indices do
                    Add(comp_colors, s.colors_list[idx]);
                od;
            fi;
            if type_char = "A" and rank = 3 and Length(comp_vecs) = 3 and Length(comp_colors) = 3 then
                a3_info := InferRealForm_A3_ByPosition(comp_vecs, comp_colors);
                if a3_info <> fail then
                    info := a3_info;
                fi;
            fi;
            if type_char = "C" and rank = 3 and Length(comp_vecs) = 3 and Length(comp_colors) = 3 then
                c3_info := InferRealForm_C3_ByPosition(comp_vecs, comp_colors);
                if c3_info <> fail then
                    info := c3_info;
                fi;
            fi;
            if type_char = "C" and rank = 2 and Length(comp_vecs) = 2 and Length(comp_colors) = 2 then
                c2_info := InferRealForm_C2_ByPosition(comp_vecs, comp_colors);
                if c2_info <> fail then
                    info := c2_info;
                fi;
            fi;
            
            if info <> fail then
                if not composite_single_like then
                    Print("        >>> 模块 2 输出: 实形式推断\n");
                    Print("        推断实形式: ", info.desc, "\n");
                    Print("        内部参数: Type=", info.type, ", Params=", info.params, "\n");
                    Print("        >>> 模块 3 输出: Noticed Orbits 匹配\n");
                fi;
                orbits := FindNoticedOrbits_Generic(info.type, info.params);
                
                if Length(orbits) > 0 then
                    if not composite_single_like then
                        Print("        发现 ", Length(orbits), " 个 Noticed Orbits\n");
                    fi;
                    Add(comp_orbits_list, orbits);
                    Add(comp_infos_list, info);
                else
                    if not composite_single_like then
                        Print("        未发现 Noticed Orbits\n");
                    fi;
                    all_comps_valid := false;
                    break;
                fi;
            else
                if not composite_single_like then
                    Print("        错误: 无法推断实形式\n");
                fi;
                all_comps_valid := false;
                break;
            fi;
        od;
        
        if all_comps_valid then
            if composite_single_like then
                Print("    >>> 模块 2 输出: 实形式推断\n");
                Print("        组件实形式: ", List(comp_infos_list, x -> x.desc), "\n");
                Print("        组件内部参数: ", List(comp_infos_list, x -> rec(Type := x.type, Params := x.params)), "\n");
            fi;
            subsystem_root_vecs := GetSubsystemRootVectorsF4(s, subsystem_roots);
            filter_info := RegisterDeltaJWKOrbitRepresentativeF4(s.idx, s.type, subsystem_root_vecs);
            Print("    >>> 模块 2.5 输出: 按 Δ_J 的 W_k-共轭类取代表\n");
            if filter_info.is_new then
                Print("        保留该 W_k·Δ_J 的代表，进入模块 3\n");
            else
                Print(Concatenation("        该 Δ_J 与索引 ", String(filter_info.first_index), " 同属一个 W_k-共轭类，跳过\n"));
            fi;
        fi;

        if all_comps_valid and filter_info.is_new then
            if composite_single_like then
                Print("    >>> 模块 3 输出: 实形式直传模块4\n");
                Print("        使用模块2组件实形式，直接进入模块4\n");
            else
                Print("    >>> 模块 3 输出: 组合轨道分析\n");
            fi;
            if IsBoundGlobal("GlobalAmbientType") and ValueGlobal("GlobalAmbientType") = "F4" then
                if not composite_single_like then
                    Print("    >>> 组合分析: F4 复合情形按组件 theta 候选统一直译，不预先枚举模块3复合轨道笛卡尔积\n");
                fi;
                Print("    ------------------------------------------\n");
                
                ordered_roots := [];
                full_ab_str := "";
                full_partition := [];
                component_payload := [];
                
                for i in [1..Length(comp_roots_list)] do
                    comp_roots := comp_roots_list[i];
                    info := comp_infos_list[i];
                    orbit := comp_orbits_list[i][1];
                    comp_colors := [];
                    for idx in s.components[i].indices do
                        Add(comp_colors, s.colors_list[idx]);
                    od;
                    
                    Append(ordered_roots, comp_roots);
                    if IsBound(orbit.partition) then
                        Append(full_partition, orbit.partition);
                    fi;
                    if IsBound(orbit.ab_diagram) then
                        for line in orbit.ab_diagram do
                            if full_ab_str <> "" then full_ab_str := Concatenation(full_ab_str, "\n"); fi;
                            full_ab_str := Concatenation(full_ab_str, line);
                        od;
                    fi;
                    
                    Add(component_payload, rec(
                        real_form_type := info.type,
                        params := info.params,
                        roots := comp_roots,
                        colors := comp_colors,
                        ab_diagram := orbit.ab_diagram
                    ));
                od;
                
                if not composite_single_like then
                    Print("    复合轨道信息:\n");
                    Print("      总划分(展示用，取各组件首个 noticed orbit): ", full_partition, "\n");
                    Print("      总 ab-diagram(展示用，取各组件首个 noticed orbit): \n", full_ab_str, "\n");
                fi;
                
                total_rank := Length(ordered_roots);
                Print("    [模块4重载] 重新加载 F4 专用 SL2 Builder\n");
                Read("../05_SL2_Triple_Construction/SL2_Triple_Builder.g");
                if IsBoundGlobal("PrintSL2TripleUnified") then
                    GlobalCompositeComponentData := component_payload;
                    ValueGlobal("PrintSL2TripleUnified")(full_ab_str, "Composite", total_rank, ordered_roots, s.colors_list);
                    Unbind(GlobalCompositeComponentData);
                else
                    PrintSL2Triple(full_ab_str, "Composite", total_rank, ordered_roots);
                fi;
                Print("\n    ------------------------------------------\n");
            else
                # 2. 组合所有组件的 Orbits (笛卡尔积)
                combined_orbits := CartesianProductList(comp_orbits_list);
                Print("    >>> 组合分析: 共生成 ", Length(combined_orbits), " 个复合轨道配置\n");
                Print("    ------------------------------------------\n");
                
                for comb in combined_orbits do
                    # comb 为 [orbit1, orbit2, ...]
                    # 将各组件的根、partition、ab-diagram 依序拼接，供 Builder 使用。
                    # 注意：保持行顺序与组件顺序一致，以确保根的消耗顺序匹配。
                    
                    # 采用策略：重排根列表 ordered_roots，使其与组件顺序一致
                    
                    ordered_roots := [];
                    full_ab_str := "";
                    full_partition := [];
                    component_payload := [];
                    
                    for i in [1..Length(comb)] do
                        orbit := comb[i];
                        comp_roots := comp_roots_list[i];
                        info := comp_infos_list[i];
                        comp_colors := [];
                        for idx in s.components[i].indices do
                            Add(comp_colors, s.colors_list[idx]);
                        od;
                        
                        # 拼接根
                        Append(ordered_roots, comp_roots);
                        
                        # 拼接 Partition
                        Append(full_partition, orbit.partition);
                        
                        # 拼接 ab-diagram
                        # orbit.ab_diagram 是字符串列表 ["aba", ...]
                        for line in orbit.ab_diagram do
                            if full_ab_str <> "" then full_ab_str := Concatenation(full_ab_str, "\n"); fi;
                            full_ab_str := Concatenation(full_ab_str, line);
                        od;
                        
                        Add(component_payload, rec(
                            real_form_type := info.type,
                            params := info.params,
                            roots := comp_roots,
                            colors := comp_colors,
                            ab_diagram := orbit.ab_diagram
                        ));
                    od;
                    
                    Print("    复合轨道信息:\n");
                    Print("      总划分: ", full_partition, "\n");
                    Print("      总 ab-diagram: \n", full_ab_str, "\n");
                    
                    # 调用 Builder (使用重排后的根列表)
                    total_rank := Length(ordered_roots);
                    # 统一入口：由构造器内部按需规范化根，再统一构造 Normal Triple
                    if IsBoundGlobal("GlobalAmbientType") and ValueGlobal("GlobalAmbientType") = "F4" then
                        Print("    [模块4重载] 重新加载 F4 专用 SL2 Builder\n");
                        Read("../05_SL2_Triple_Construction/SL2_Triple_Builder.g");
                    fi;
                    if IsBoundGlobal("PrintSL2TripleUnified") then
                        GlobalCompositeComponentData := component_payload;
                        ValueGlobal("PrintSL2TripleUnified")(full_ab_str, "Composite", total_rank, ordered_roots, s.colors_list);
                        Unbind(GlobalCompositeComponentData);
                    else
                        PrintSL2Triple(full_ab_str, "Composite", total_rank, ordered_roots);
                    fi;
                    Print("\n    ------------------------------------------\n");
                od;
            fi;
        fi;
        
    else
        # 原有的单子系统处理逻辑
        
        # 解析子系统类型 (例如 "A2" -> type="A", rank=2；"B4/C4" 取 "B4")
        raw_type := s.type;
        slash_pos := Position(raw_type, '/');
        if slash_pos <> fail and slash_pos > 1 then
            raw_type := raw_type{[1..slash_pos-1]};
        fi;
        type_char := raw_type{[1]};
        if Length(raw_type) > 1 then
            rank_str := raw_type{[2..Length(raw_type)]};
            if IsDigitChar(rank_str[1]) then
                rank := Int(rank_str);
                
                # 推断实形式
                if type_char = "B" and rank = 4 and IsBound(s.colors_list) then
                    info := InferRealForm_B4_FromColors(s.colors_list);
                else
                    info := InferRealForm(type_char, rank, s.black_nodes, s.white_nodes);
                fi;
                if type_char = "C" and rank = 2 and IsBound(s.roots_vecs) and IsBound(s.colors_list) and Length(s.roots_vecs) = 2 and Length(s.colors_list) = 2 then
                    c2_info := InferRealForm_C2_ByPosition(s.roots_vecs, s.colors_list);
                    if c2_info <> fail then
                        info := c2_info;
                    fi;
                fi;
                if type_char = "C" and rank = 3 and IsBound(s.roots_vecs) and IsBound(s.colors_list) and Length(s.roots_vecs) = 3 and Length(s.colors_list) = 3 then
                    c3_info := InferRealForm_C3_ByPosition(s.roots_vecs, s.colors_list);
                    if c3_info <> fail then
                        info := c3_info;
                    fi;
                fi;
                
                if info <> fail then
                    Print("    >>> 模块 2 输出: 实形式推断\n");
                    Print("        实形式: ", info.desc, "\n");
                    Print("        内部参数: Type=", info.type, ", Params=", info.params, "\n");
                    GlobalModule2CurrentRealFormInfo := rec(type := info.type, params := info.params);
                    Print("    >>> 模块 2.5 输出: 按 Δ_J 的 W_k-共轭类取代表\n");
                    subsystem_root_vecs := GetSubsystemRootVectorsF4(s, subsystem_roots);
                    filter_info := RegisterDeltaJWKOrbitRepresentativeF4(s.idx, s.type, subsystem_root_vecs);
                    if filter_info.is_new then
                        Print("        保留该 W_k·Δ_J 的代表，进入模块 3\n");
                        Print("    ------------------------------------------\n");
                        
                        # 调用 Noticed Orbits 分析
                        if type_char = "B" and rank = 4 and IsBound(s.colors_list) then
                            GlobalColorsListForDiagram := s.colors_list;
                        fi;
                        if IsBoundGlobal("PrintNoticedOrbits_Generic") then
                            ValueGlobal("PrintNoticedOrbits_Generic")(info.type, info.params, subsystem_roots);
                        else
                            Print("    警告: 未加载 Noticed_Orbits 模块，跳过模块 3 输出\n");
                        fi;
                        AppendPipelineLog("CALL Noticed\n");
                        if IsBound(GlobalColorsListForDiagram) then
                            Unbind(GlobalColorsListForDiagram);
                        fi;
                        if type_char = "B" and rank = 4 then
                            Print("    >>> 追加 F4 诊断输出\n");
                            AppendPipelineLog("DIAG F4\n");
                            # F4诊断输出已在PrintNoticedOrbits_Generic中完成
                        fi;
                    else
                        Print(Concatenation("        该 Δ_J 与索引 ", String(filter_info.first_index), " 同属一个 W_k-共轭类，跳过\n"));
                    fi;
                else
                    Print("    错误: 无法推断实形式\n");
                fi;
            else
                Print("    错误: 无法解析秩 (可能是复合类型但未正确识别组件)\n");
            fi;
        else
            Print("    错误: 类型字符串太短\n");
        fi;
    fi;
    if IsBound(GlobalModule2CurrentRealFormInfo) then
        Unbind(GlobalModule2CurrentRealFormInfo);
    fi;
    Print("    ==========================================\n");
od;

if IsBoundGlobal("GlobalNormalTripleTotalChecks") then
    EmitSummary("Module4 NormalTriple", Concatenation(String(GlobalNormalTriplePassChecks), "/", String(GlobalNormalTripleTotalChecks)));
fi;
if IsBoundGlobal("GlobalF4RawLabelStats") then
    EmitSummary("Module5 RawCompactDigits", GlobalF4RawLabelStats);
elif IsBoundGlobal("GlobalF4LabelStats") then
    EmitSummary("Module5 RawCompactDigits", GlobalF4LabelStats);
fi;
EmitSummaryIfGlobal("Module5 TerminalKOrbitLabels", "GlobalF4TerminalLabelStats");
if IsBoundGlobal("GlobalF4LabelStats") then
    EmitSummary("Module5 KOrbitLabels", ValueGlobal("GlobalF4LabelStats"));
    EmitSummaryIfGlobal("Module5 FallbackCount", "GlobalF4LabelFallbackCount");
fi;
EmitSummaryIfGlobal("Module5 MappingTotal", "GlobalF4MappingTotal");
EmitSummaryIfGlobal("Module5 AllowedViolationCount", "GlobalF4AllowedViolationCount");
EmitSummaryIfGlobal("Module5 TerminalConvergedCount", "GlobalF4TerminalConvergedCount");
EmitSummaryIfGlobal("Module5 TerminalViolationCount", "GlobalF4TerminalViolationCount");
EmitSummaryIfGlobal("Module5 IdempotentConvergedCount", "GlobalF4ConvergenceIdempotentCount");
EmitSummaryIfGlobal("Module5 IdempotentViolationCount", "GlobalF4ConvergenceIdempotentViolationCount");
EmitSummaryIfGlobal("Module5 AffineConstraintTotal", "GlobalF4AffineConstraintTotal");
EmitSummaryIfGlobal("Module5 AffineConstraintViolationCount", "GlobalF4AffineConstraintViolationCount");
EmitSummaryIfGlobal("Module5 CompactWeightSetTotal", "GlobalF4CompactWeightSetTotal");
EmitSummaryIfGlobal("Module5 CompactWeightSetViolationCount", "GlobalF4CompactWeightSetViolationCount");
if IsBoundGlobal("F4RunConvergenceStressTest") then
    stress_res := ValueGlobal("F4RunConvergenceStressTest")();
    EmitSummary("Module5 StressTotal", stress_res.total);
    EmitSummary("Module5 StressAllowedViolationCount", stress_res.allowed_bad);
    EmitSummary("Module5 StressTerminalViolationCount", stress_res.terminal_bad);
    EmitSummary("Module5 StressIdempotentViolationCount", stress_res.idempotent_bad);
fi;
assert_m4_all_pass := GetGlobalOrDefault("GlobalRequireModule4AllPass", false);
assert_m5_labels_allowed := GetGlobalOrDefault("GlobalRequireModule5LabelsInAllowed", true);
assert_m5_no_fallback := GetGlobalOrDefault("GlobalRequireModule5NoFallback", false);
assert_m5_terminal_only := GetGlobalOrDefault("GlobalRequireModule5TerminalOnly", true);
assert_m5_idempotent := GetGlobalOrDefault("GlobalRequireModule5Idempotent", true);
assert_m5_stress_pass := GetGlobalOrDefault("GlobalRequireModule5StressPass", true);
assert_m5_affine_constraint := GetGlobalOrDefault("GlobalRequireModule5AffineConstraint", true);
assert_m5_compact_weight_set := GetGlobalOrDefault("GlobalRequireModule5CompactWeightSet", true);
assert_failures := [];
if assert_m4_all_pass and IsBoundGlobal("GlobalNormalTripleTotalChecks") then
    if GlobalNormalTriplePassChecks <> GlobalNormalTripleTotalChecks then
        Add(assert_failures, Concatenation("模块4断言失败: NormalTriple通过数=", String(GlobalNormalTriplePassChecks), "/", String(GlobalNormalTripleTotalChecks)));
    fi;
fi;
if assert_m5_labels_allowed and IsBoundGlobal("F4TerminalKOrbitLabels") then
    allowed_labels := ValueGlobal("F4TerminalKOrbitLabels")();
    if IsBoundGlobal("GlobalF4LabelStats") then
        for pair in GlobalF4LabelStats do
            if Position(allowed_labels, pair[1]) = fail then
                Add(assert_failures, Concatenation("模块5断言失败: 出现目标列表外标签 ", String(pair[1])));
            fi;
        od;
    fi;
fi;
if assert_m5_no_fallback and IsBoundGlobal("GlobalF4LabelFallbackCount") then
    if GlobalF4LabelFallbackCount <> 0 then
        Add(assert_failures, Concatenation("模块5断言失败: FallbackCount=", String(GlobalF4LabelFallbackCount), " (期望0)"));
    fi;
fi;
if assert_m5_terminal_only and IsBoundGlobal("GlobalF4TerminalLabelStats") then
    allowed_labels := [];
    if IsBoundGlobal("F4TerminalKOrbitLabels") then
        allowed_labels := ValueGlobal("F4TerminalKOrbitLabels")();
    fi;
    for pair in GlobalF4TerminalLabelStats do
        if Position(allowed_labels, pair[1]) = fail then
            Add(assert_failures, Concatenation("模块5断言失败: 非终态标签 ", String(pair[1])));
        fi;
    od;
fi;
if assert_m5_terminal_only and IsBoundGlobal("GlobalF4TerminalViolationCount") then
    if GlobalF4TerminalViolationCount <> 0 then
        Add(assert_failures, Concatenation("模块5断言失败: TerminalViolationCount=", String(GlobalF4TerminalViolationCount), " (期望0)"));
    fi;
fi;
if assert_m5_idempotent and IsBoundGlobal("GlobalF4ConvergenceIdempotentViolationCount") then
    if GlobalF4ConvergenceIdempotentViolationCount <> 0 then
        Add(assert_failures, Concatenation("模块5断言失败: IdempotentViolationCount=", String(GlobalF4ConvergenceIdempotentViolationCount), " (期望0)"));
    fi;
fi;
if assert_m5_stress_pass and IsBoundGlobal("GlobalF4StressAllowedViolationCount") then
    if GlobalF4StressAllowedViolationCount <> 0 then
        Add(assert_failures, Concatenation("模块5断言失败: StressAllowedViolationCount=", String(GlobalF4StressAllowedViolationCount), " (期望0)"));
    fi;
    if GlobalF4StressTerminalViolationCount <> 0 then
        Add(assert_failures, Concatenation("模块5断言失败: StressTerminalViolationCount=", String(GlobalF4StressTerminalViolationCount), " (期望0)"));
    fi;
    if GlobalF4StressIdempotentViolationCount <> 0 then
        Add(assert_failures, Concatenation("模块5断言失败: StressIdempotentViolationCount=", String(GlobalF4StressIdempotentViolationCount), " (期望0)"));
    fi;
fi;
if assert_m5_affine_constraint and IsBoundGlobal("GlobalF4AffineConstraintViolationCount") then
    if GlobalF4AffineConstraintViolationCount <> 0 then
        Add(assert_failures, Concatenation("模块5断言失败: AffineConstraintViolationCount=", String(GlobalF4AffineConstraintViolationCount), " (期望0)"));
    fi;
fi;
if assert_m5_compact_weight_set and IsBoundGlobal("GlobalF4CompactWeightSetViolationCount") then
    if GlobalF4CompactWeightSetViolationCount <> 0 then
        Add(assert_failures, Concatenation("模块5断言失败: CompactWeightSetViolationCount=", String(GlobalF4CompactWeightSetViolationCount), " (期望0)"));
    fi;
fi;
if Length(assert_failures) = 0 then
    Print("ASSERT SUMMARY: PASS\n");
    AppendPipelineLog("ASSERT SUMMARY PASS\n");
else
    Print("ASSERT SUMMARY: FAIL\n");
    AppendPipelineLog("ASSERT SUMMARY FAIL\n");
    for err in assert_failures do
        Print("  - ", err, "\n");
        AppendPipelineLog(Concatenation("  - ", err, "\n"));
    od;
    Print("\n全流程分析完成。\n");
    QUIT_GAP(1);
fi;
Print("\n全流程分析完成。\n");
quit;
