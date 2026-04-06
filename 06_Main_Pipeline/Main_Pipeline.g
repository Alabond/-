#############################################################################
# 文件: Main_Pipeline.g
# 描述: 通用主流程执行器。
# 作用: 读取各模块、枚举实验、筛选子系统、推断实形式、匹配 Noticed Orbit、
#       调用模块 4 构造 normal triple，并最终汇总断言结果。
# 说明: 该文件本身不固定某个具体实验配置，通常由 G2_Pipeline.g 之类的入口
#       先设置全局变量，再通过 Read("Main_Pipeline.g") 调用。
# 输出: 当前版本只保留一个完整结果文件，不再单独维护摘要日志文件。
#############################################################################

# 启动标记，方便从日志中观察主流程是否开始执行。
Print("BOOT\n");

# 允许外部覆盖输出文件名
if not IsBound(GlobalOutputDir) then
    GlobalOutputDir := "outputs";
fi;
out_file := Concatenation(GlobalOutputDir, "/g2(2)实验运行结果.txt");
if IsBound(GlobalOutFile) then
    out_file := GlobalOutFile;
fi;
LogTo(out_file);

if not IsBound(GlobalAmbientType) then
    GlobalAmbientType := "G2";
fi;

# 首先加载基础工具模块
Print("LOAD 02_Subsystem_Classification\n");
Read("../02_Subsystem_Classification/Subsystem_Classifier.g");
Print("LOAD 04_Orbit_Analysis/K_Orbit_Classifier\n");
Read("../04_Orbit_Analysis/K_Orbit_Classifier.g");

# 然后加载其他模块
Print("LOAD 01_Subsystem_Generation\n");
Read("../01_Subsystem_Generation/G2_Extended_Subsystem_Lib.g");
Print("LOAD 04_Orbit_Analysis/Noticed_Orbits\n");
Read("../04_Orbit_Analysis/Noticed_Orbits.g");
Print("LOAD 03_Real_Form_Inference\n");
Read("../03_Real_Form_Inference/Real_Form_Inference.g");
Print("LOAD 05_SL2_Triple_Construction\n");
Read("../05_SL2_Triple_Construction/SL2_Triple_Builder.g");

Print("START\n");

# 输出一条统一格式的摘要信息。
# 用于汇总 normal triple 通过率、断言结果等总览指标。
EmitSummary := function(key, value)
    Print("SUMMARY ", key, ": ", value, "\n");
end;;

# 读取一个全局变量；若未定义则返回给定默认值。
# 用于让配置层按需覆盖主流程中的默认参数。
GetGlobalOrDefault := function(global_name, default_value)
    if IsBoundGlobal(global_name) then
        return ValueGlobal(global_name);
    fi;
    return default_value;
end;;

GlobalNormalTripleTotalChecks := 0;
GlobalNormalTriplePassChecks := 0;

# 2. 步骤 1: G2 子系统筛选与配置
Print("\n==================================================\n");
Print("步骤 1: 子系统选取与配置\n");
Print("==================================================\n");
# 允许外部覆盖简介文本
if IsBound(GlobalIntroText) then
    Print(GlobalIntroText, "\n");
else
    Print("研究对象: G2(2) 两组实验 -> A2(一黑一白), A1+A1(两黑)\n");
fi;

# 用户配置区域（允许外部注入实验；默认仅运行 G2 的两个实验）
if IsBound(GlobalExperiments) then
    experiments := GlobalExperiments;
else
    experiments := [
        rec(
            desc := "G2 实验 1: 长-长 一黑一白 -> A2",
            roots := rec(
                mode := "SCAN",
                template := [[0,1], [0,1]],
                target_black := 1,
                target_white := 1
            ),
            colors := fail
        ),
        rec(
            desc := "G2 实验 2: 两黑 正交 -> A1 + A1",
            roots := rec(
                mode := "SCAN",
                template := [[1,0], [0,1]],
                target_black := 2,
                target_white := 0
            ),
            colors := fail
        )
    ];
fi;

all_subsystems := [];

# 逐个实验配置调用模块 1 进行子系统扫描或直接分类。
# 每个实验可能返回多个候选子系统，最后统一汇总到 all_subsystems。
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
        elif IsBound(exp.colors) and exp.colors <> fail then
            Print("  Target colors: ", exp.colors, "\n");
        fi;
    elif exp.roots <> fail then
        Print("  Roots: ", exp.roots, "\n");
        Print("  Colors: ", exp.colors, "\n");
    else
        Print("  Mode: Auto Scan\n");
    fi;
    
    subs := FindG2Subsystems(exp.roots, exp.colors);
    Print("  发现 ", Length(subs), " 个子系统配置。\n");
    
    # 为子系统添加实验描述标记，以便区分
    for s in subs do
        s.exp_desc := exp.desc;
    od;
    
    Append(all_subsystems, subs);
od;

subsystems := all_subsystems;

# 文本对齐辅助（右侧填充空格至固定宽度）
# 把字符串右侧补空格到固定宽度。
# 用于控制命令行表格输出的列对齐效果。
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
    
    Print(
        PadRight(String(s.idx), 6), "| ",
        PadRight(roots_display, 34), "| ",
        PadRight(Concatenation(color_desc, black_pos_desc), 20), "| ",
        PadRight(s.type, 10), "\n" # 对应增加宽度
    );
od;

# 3. 步骤 2: 自动分析每个唯一子系统类型的 Noticed Orbits
Print("\n==================================================\n");
Print("步骤 2: 自动分析 Noticed Orbits\n");
Print("==================================================\n");

# 辅助函数: 计算列表族的笛卡尔积
# 例如 [[a,b],[1,2]] -> [[a,1],[a,2],[b,1],[b,2]]
# 用于复合子系统场景下，把各组件的 noticed orbit 列表做全组合。
# 返回的每个元素都是“每个组件各选一个 orbit”组成的组合配置。
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

CanonicalDeltaJWKOrbitKeyG2 := function(type_str, roots_list)
    local generators, words, gcj, parsed_roots, best_key, word, entries, i, v, k, entry, key, sep;
    generators := [[1, 0], [-3, -2]];
    words := [[], [1], [2], [1, 2]];
    gcj := fail;
    if IsBoundGlobal("BuildGCJFromJ_G2") then
        gcj := ValueGlobal("BuildGCJFromJ_G2")(roots_list);
    fi;
    if gcj <> fail and IsRecord(gcj) and IsBound(gcj.delta_roots) then
        parsed_roots := List(gcj.delta_roots, r -> ValueGlobal("ParseRootString")(r));
    else
        parsed_roots := List(roots_list, r -> ValueGlobal("ParseRootString")(r));
    fi;
    best_key := fail;
    for word in words do
        entries := [];
        for i in [1..Length(parsed_roots)] do
            v := [parsed_roots[i][1], parsed_roots[i][2]];
            for k in word do
                v := v - (2 * ValueGlobal("InnerProduct")(v, generators[k]) / ValueGlobal("InnerProduct")(generators[k], generators[k])) * generators[k];
            od;
            entry := ValueGlobal("FormatRoot")(v);
            Add(entries, entry);
        od;
        entries := SortedList(entries);
        key := Concatenation(type_str, "|", String(Length(entries)), "|");
        sep := "";
        for entry in entries do
            key := Concatenation(key, sep, entry);
            sep := ";";
        od;
        if best_key = fail or key < best_key then
            best_key := key;
        fi;
    od;
    return best_key;
end;;

RegisterDeltaJWKOrbitRepresentativeG2 := function(subsystem_index, type_str, roots_list)
    local key, item;
    if not IsBoundGlobal("GlobalDeltaJWKOrbitSeen") then
        GlobalDeltaJWKOrbitSeen := [];
    fi;
    key := CanonicalDeltaJWKOrbitKeyG2(type_str, roots_list);
    for item in GlobalDeltaJWKOrbitSeen do
        if item.key = key then
            return rec(is_new := false, first_index := item.first_index, key := key);
        fi;
    od;
    Add(GlobalDeltaJWKOrbitSeen, rec(key := key, first_index := subsystem_index));
    return rec(is_new := true, first_index := subsystem_index, key := key);
end;;

# 遍历每个候选子系统并推进完整分析链。
# 若子系统是复合型，则逐组件分析后再做笛卡尔积组合；
# 若是单一型，则直接走“实形式推断 -> noticed orbit -> sl2-triple”主链。
GlobalDeltaJWKOrbitSeen := [];
GlobalDeltaJWKOrbitExperimentTag := fail;
for s in subsystems do
    if IsBound(s.exp_desc) and GlobalDeltaJWKOrbitExperimentTag <> s.exp_desc then
        GlobalDeltaJWKOrbitSeen := [];
        GlobalDeltaJWKOrbitExperimentTag := s.exp_desc;
    fi;
    Print("\n>>> 模块 1 输出: 子系统分析 (索引 ", s.idx, ")\n");
    if IsBound(s.exp_desc) then
        Print("    实验: ", s.exp_desc, "\n");
    fi;
    Print("    类型: ", s.type, "\n");
    
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
        Print("    >>> 检测到复合子系统，包含 ", Length(s.components), " 个组件\n");
        
        # 1. 对每个组件分别进行实形式推断和 Noticed Orbit 分析
        comp_roots_list := [];
        comp_infos_list := [];
        
        all_comps_valid := true;
        
        for comp in s.components do
            # comp.indices: 组件在 subsystem_roots 中的索引列表
            # 提取该组件的根
            comp_roots := [];
            for idx in comp.indices do
                Add(comp_roots, subsystem_roots[idx]);
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
            type_char := comp.type{[1]};
            rank_str := comp.type{[2..Length(comp.type)]};
            rank := Int(rank_str);
            
            Print("    --- 组件: ", comp.type, " (根: ", comp_roots, ", 颜色: ", comp_black, "B ", comp_white, "W) ---\n");
            
            # 推断实形式
            info := InferRealForm(type_char, rank, comp_black, comp_white);
            
            if info <> fail then
                Print("        >>> 模块 2 输出: 实形式推断\n");
                Print("        推断实形式: ", info.desc, "\n");
                Print("        内部参数: Type=", info.type, ", Params=", info.params, "\n");
                Add(comp_infos_list, info);
            else
                Print("        错误: 无法推断实形式\n");
                all_comps_valid := false;
                break;
            fi;
        od;
        
        if all_comps_valid then
            filter_info := RegisterDeltaJWKOrbitRepresentativeG2(s.idx, s.type, subsystem_roots);
            Print("    >>> 模块 2.5 输出: 按 Δ_J 的 W_k-共轭类取代表\n");
            if filter_info.is_new then
                Print("        保留该 W_k·Δ_J 的代表，进入模块 3\n");
            else
                Print(Concatenation("        该 Δ_J 与索引 ", String(filter_info.first_index), " 同属一个 W_k-共轭类，跳过\n"));
            fi;
        fi;

        if all_comps_valid and filter_info.is_new then
            comp_orbits_list := [];
            for i in [1..Length(comp_infos_list)] do
                info := comp_infos_list[i];
                Print("        >>> 模块 3 输出: Noticed Orbits 匹配\n");
                orbits := FindNoticedOrbits_Generic(info.type, info.params);
                if Length(orbits) > 0 then
                    Print("        发现 ", Length(orbits), " 个 Noticed Orbits\n");
                    Add(comp_orbits_list, orbits);
                else
                    Print("        未发现 Noticed Orbits\n");
                    all_comps_valid := false;
                    break;
                fi;
            od;
        fi;

        if all_comps_valid and filter_info.is_new then
            Print("    >>> 模块 3 输出: 组合轨道分析\n");
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
                if IsBoundGlobal("PrintSL2TripleUnified") then
                    GlobalCompositeComponentData := component_payload;
                    if Length(component_payload) = 1 then
                        GlobalKOrbitSubsystemTypeOverride := s.type;
                        GlobalModule4CurrentRealFormInfo := rec(type := component_payload[1].real_form_type, params := component_payload[1].params);
                    fi;
                    ValueGlobal("PrintSL2TripleUnified")(full_ab_str, "Composite", total_rank, ordered_roots, s.colors_list);
                    Unbind(GlobalCompositeComponentData);
                    if IsBoundGlobal("GlobalKOrbitSubsystemTypeOverride") then
                        Unbind(GlobalKOrbitSubsystemTypeOverride);
                    fi;
                    if IsBoundGlobal("GlobalModule4CurrentRealFormInfo") then
                        Unbind(GlobalModule4CurrentRealFormInfo);
                    fi;
                else
                    PrintSL2Triple(full_ab_str, "Composite", total_rank, ordered_roots);
                fi;
                Print("\n    ------------------------------------------\n");
            od;
        fi;
        
    else
        # 原有的单子系统处理逻辑
        
        # 解析子系统类型 (例如 "A2" -> type="A", rank=2)
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
                info := InferRealForm(type_char, rank, s.black_nodes, s.white_nodes);
                
                if info <> fail then
                    Print("    >>> 模块 2 输出: 实形式推断\n");
                    Print("        实形式: ", info.desc, "\n");
                    Print("        内部参数: Type=", info.type, ", Params=", info.params, "\n");
                    Print("    >>> 模块 2.5 输出: 按 Δ_J 的 W_k-共轭类取代表\n");
                    filter_info := RegisterDeltaJWKOrbitRepresentativeG2(s.idx, s.type, subsystem_roots);
                    if filter_info.is_new then
                        Print("        保留该 W_k·Δ_J 的代表，进入模块 3\n");
                        Print("    ------------------------------------------\n");
                        if IsBoundGlobal("PrintNoticedOrbits_Generic") then
                            ValueGlobal("PrintNoticedOrbits_Generic")(info.type, info.params, subsystem_roots);
                        else
                            Print("    警告: 未加载 Noticed_Orbits 模块，跳过模块 3 输出\n");
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
    Print("    ==========================================\n");
od;

# 汇总模块 4 的 normal triple 校验统计，并在需要时执行最终断言。
# 若断言失败则以非零状态退出 GAP，方便外部脚本识别失败。
if IsBoundGlobal("GlobalNormalTripleTotalChecks") then
    EmitSummary("Module4 NormalTriple", Concatenation(String(GlobalNormalTriplePassChecks), "/", String(GlobalNormalTripleTotalChecks)));
fi;
assert_m4_all_pass := GetGlobalOrDefault("GlobalRequireModule4AllPass", true);
assert_failures := [];
if assert_m4_all_pass and IsBoundGlobal("GlobalNormalTripleTotalChecks") then
    if GlobalNormalTriplePassChecks <> GlobalNormalTripleTotalChecks then
        Add(assert_failures, Concatenation("模块4断言失败: NormalTriple通过数=", String(GlobalNormalTriplePassChecks), "/", String(GlobalNormalTripleTotalChecks)));
    fi;
fi;
if Length(assert_failures) = 0 then
    Print("ASSERT SUMMARY: PASS\n");
else
    Print("ASSERT SUMMARY: FAIL\n");
    for err in assert_failures do
        Print("  - ", err, "\n");
    od;
    Print("\n全流程分析完成。\n");
    QUIT_GAP(1);
fi;
Print("\n全流程分析完成。\n");
quit;
