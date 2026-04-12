#############################################################################
# 文件: Noticed_Orbits.g
# 描述: 经典实李代数 Noticed Nilpotent Orbits 与 ab-diagram 生成
# 参考: Noel, "Nilpotent orbits and theta-stable parabolic subalgebras" (1994)
#############################################################################

# =============================================================================
# 第一部分: Noticed Orbit 判定逻辑
# =============================================================================

GlobalModule4CurrentRealFormInfo := fail;

# 1. Type A: sl(n, R) (Split Form)
# -----------------------------------------------------------------------------
# Noel Thm 4.2.1: Partitions of n with distinct parts
IsNoticed_TypeA_Split := function(partition)
    if Length(Set(partition)) < Length(partition) then return false; fi;
    return true;
end;;

# 1b. Type A: su(p, q) (Non-Split Form)
# -----------------------------------------------------------------------------
# Noel Thm 4.2.1:
# - su(p, p+1): Exactly 1 noticed orbit (one-row signed Young diagram of size 2p+1)
# - su(p, p): Exactly 2 noticed orbits (one-row signed Young diagram of size 2p)
# - |p-q| >= 2: No noticed orbits
IsNoticed_TypeA_SU_PQ := function(partition, p, q)
    local n;
    n := p + q;
    # 必须是单行 (partition = [n])
    if Length(partition) <> 1 or partition[1] <> n then
        return false;
    fi;
    # 检查 p, q 条件
    if AbsInt(p - q) <= 1 then return true; else return false; fi;
end;;

# 1c. Type A: su*(2n)
# -----------------------------------------------------------------------------
# Noel Thm 4.2.1: No non-zero noticed orbits
IsNoticed_TypeA_SU_Star := function(partition)
    return false;
end;;

# 2. Type C: sp(n, R)
# -----------------------------------------------------------------------------
# Noel Thm 4.2.3:
# - All rows even
# - Each row length can repeat at most twice
IsNoticed_TypeC_Split := function(partition)
    local len, counts, c;
    for len in partition do
        if len mod 2 <> 0 then return false; fi;
    od;
    counts := Collected(partition);
    for c in counts do
        if c[2] > 2 then return false; fi;
    od;
    return true;
end;;

# 2b. Type C: sp(p, q)
# -----------------------------------------------------------------------------
# Noel Thm 4.2.3: No non-zero noticed orbits
IsNoticed_TypeC_SP_PQ := function(partition)
    return false;
end;;

# 3. Type B/D: so(p, q)
# -----------------------------------------------------------------------------
# Noel Thm 4.2.2:
# - All rows odd
# - Each row length can repeat at most twice
IsNoticed_TypeBD_Orthogonal := function(partition)
    local len, counts, c;
    for len in partition do
        if len mod 2 = 0 then return false; fi; # 必须全为奇数
    od;
    counts := Collected(partition);
    for c in counts do
        if c[2] > 2 then return false; fi;
    od;
    return true;
end;;

# 3b. Type D: so*(2n)
# -----------------------------------------------------------------------------
# Noel Thm 4.2.2: No non-zero noticed orbits
IsNoticed_TypeD_SO_Star := function(partition)
    return false;
end;;

# ab-diagram 生成 (通用)
GenerateABDiagram_TypeA := function(partition)
    local diagram, row, len, j, s;
    diagram := [];
    for len in partition do
        s := "";
        for j in [1..len] do
            if j mod 2 = 1 then s := Concatenation(s, "a");
            else s := Concatenation(s, "b"); fi;
        od;
        Add(diagram, s);
    od;
    return diagram;
end;;

GenerateABDiagram_TypeA_StartB := function(partition)
    local diagram, len, j, s;
    diagram := [];
    for len in partition do
        s := "";
        for j in [1..len] do
            if j mod 2 = 1 then s := Concatenation(s, "b");
            else s := Concatenation(s, "a"); fi;
        od;
        Add(diagram, s);
    od;
    return diagram;
end;;

# =============================================================================
# 第二部分: 主查询接口
# =============================================================================

# 通用查询入口
# real_form_type: 
#   "sl_R" (Split Type A)
#   "su_pq" (Type A)
#   "su_star" (Type A)
#   "sp_R" (Split Type C)
#   "sp_pq" (Type C)
#   "so_pq" (Type B/D Split & Non-Split)
#   "so_star" (Type D)
#
# params: [rank] or [p, q] or [n]
FindNoticedOrbits_Generic := function(real_form_type, params)
    local N, parts, p, result, ab, r_p, r_q, ab_alt;
    
    result := [];
    
    if real_form_type = "sl_R" then
        # params: [rank] -> sl(rank+1, R)
        N := params[1] + 1;
        parts := Partitions(N);
        for p in parts do
            if IsNoticed_TypeA_Split(p) then
                ab := GenerateABDiagram_TypeA(p);
                Add(result, rec(partition := p, ab_diagram := ab, desc := "sl(n, R)"));
            fi;
        od;
        
    elif real_form_type = "su_pq" then
        # params: [p, q]
        r_p := params[1];
        r_q := params[2];
        N := r_p + r_q;
        parts := Partitions(N);
        for p in parts do
            if IsNoticed_TypeA_SU_PQ(p, r_p, r_q) then
                ab := GenerateABDiagram_TypeA(p);
                if r_p = r_q and Length(p) = 1 and p[1] = N then
                    Add(result, rec(partition := p, ab_diagram := ab, desc := Concatenation("su(", String(r_p), ",", String(r_q), ")")));
                    ab_alt := GenerateABDiagram_TypeA_StartB(p);
                    Add(result, rec(partition := p, ab_diagram := ab_alt, desc := Concatenation("su(", String(r_p), ",", String(r_q), ")")));
                else
                    Add(result, rec(partition := p, ab_diagram := ab, desc := Concatenation("su(", String(r_p), ",", String(r_q), ")")));
                fi;
            fi;
        od;
        
    elif real_form_type = "sp_R" then
        # params: [rank] -> sp(2*rank, R)
        N := 2 * params[1];
        parts := Partitions(N);
        for p in parts do
            if IsNoticed_TypeC_Split(p) then
                ab := GenerateABDiagram_TypeA(p);
                Add(result, rec(partition := p, ab_diagram := ab, desc := "sp(n, R)"));
            fi;
        od;
        
    elif real_form_type = "so_pq" then
        # params: [p, q]
        r_p := params[1];
        r_q := params[2];
        N := r_p + r_q;
        parts := Partitions(N);
        for p in parts do
            if IsNoticed_TypeBD_Orthogonal(p) then
                ab := GenerateABDiagram_TypeA(p);
                Add(result, rec(partition := p, ab_diagram := ab, desc := Concatenation("so(", String(r_p), ",", String(r_q), ")")));
            fi;
        od;
    
    # 其他无 Noticed Orbits 的情况
    elif real_form_type = "su_star" or real_form_type = "sp_pq" or real_form_type = "so_star" then
        return [];
    fi;
    
    return result;
end;;

# 修改接口调用
PrintNoticedOrbits_Generic := function(real_form_type, params, subsystem_roots)
    local orbits, o, ab_str, i, total_edges, row, expanded_roots, n_base;
    if IsBoundGlobal("GlobalAmbientType") and ValueGlobal("GlobalAmbientType") = "F4" and subsystem_roots <> fail then
        Print("    >>> 模块 3 输出: 实形式直传模块4\n");
        if IsBoundGlobal("AppendRetainedMirrorText") then
            ValueGlobal("AppendRetainedMirrorText")("    >>> 模块 3 输出: 实形式直传模块4\n");
        fi;
        Print("        使用模块2实形式: RealForm=", real_form_type, ", Params=", params, "\n");
        if IsBoundGlobal("AppendRetainedMirrorText") then
            ValueGlobal("AppendRetainedMirrorText")(Concatenation("        使用模块2实形式: RealForm=", real_form_type, ", Params=", String(params), "\n"));
        fi;
        if IsBoundGlobal("GlobalModule4CurrentRealFormInfo") then
            Unbind(GlobalModule4CurrentRealFormInfo);
        fi;
        GlobalModule4CurrentRealFormInfo := rec(type := real_form_type, params := params);
        Print("        直接进入模块4，由 theta-orbit 路径枚举 noticed 轨道\n");
        if IsBoundGlobal("AppendRetainedMirrorText") then
            ValueGlobal("AppendRetainedMirrorText")("        直接进入模块4，由 theta-orbit 路径枚举 noticed 轨道\n");
        fi;
        if IsBoundGlobal("GlobalAmbientType") and ValueGlobal("GlobalAmbientType") = "F4" then
            Print("        [模块4重载] 重新加载 F4 专用 SL2 Builder\n");
            if IsBoundGlobal("AppendRetainedMirrorText") then
                ValueGlobal("AppendRetainedMirrorText")("        [模块4重载] 重新加载 F4 专用 SL2 Builder\n");
            fi;
            Read("../05_SL2_Triple_Construction/SL2_Triple_Builder.g");
        fi;
        if IsBoundGlobal("PrintSL2TripleUnified") then
            ValueGlobal("PrintSL2TripleUnified")("", real_form_type, params[1], subsystem_roots, fail);
        elif IsBoundGlobal("PrintSL2Triple") then
            ValueGlobal("PrintSL2Triple")("", real_form_type, params[1], subsystem_roots, fail);
        fi;
        if IsBoundGlobal("GlobalModule4CurrentRealFormInfo") then
            Unbind(GlobalModule4CurrentRealFormInfo);
        fi;
        return;
    fi;
    orbits := FindNoticedOrbits_Generic(real_form_type, params);
    
    if Length(orbits) = 0 then
        Print("    >>> 模块 3 输出: Noticed Orbit 分析\n");
        if IsBoundGlobal("AppendRetainedMirrorText") then
            ValueGlobal("AppendRetainedMirrorText")("    >>> 模块 3 输出: Noticed Orbit 分析\n");
        fi;
        Print("        未发现 Noticed Nilpotent 轨道。\n");
        if IsBoundGlobal("AppendRetainedMirrorText") then
            ValueGlobal("AppendRetainedMirrorText")("        未发现 Noticed Nilpotent 轨道。\n");
        fi;
    else
        Print("    >>> 模块 3 输出: Noticed Orbit 分析\n");
        if IsBoundGlobal("AppendRetainedMirrorText") then
            ValueGlobal("AppendRetainedMirrorText")("    >>> 模块 3 输出: Noticed Orbit 分析\n");
        fi;
        Print("        查询参数: RealForm=", real_form_type, ", Params=", params, "\n");
        if IsBoundGlobal("AppendRetainedMirrorText") then
            ValueGlobal("AppendRetainedMirrorText")(Concatenation("        查询参数: RealForm=", real_form_type, ", Params=", String(params), "\n"));
        fi;
        Print("        共发现 ", Length(orbits), " 个 Noticed Nilpotent 轨道。\n");
        if IsBoundGlobal("AppendRetainedMirrorText") then
            ValueGlobal("AppendRetainedMirrorText")(Concatenation("        共发现 ", String(Length(orbits)), " 个 Noticed Nilpotent 轨道。\n"));
        fi;
        
        for o in orbits do
            Print("        --------------------------------------\n");
            if IsBoundGlobal("AppendRetainedMirrorText") then
                ValueGlobal("AppendRetainedMirrorText")("        --------------------------------------\n");
            fi;
            Print("        轨道信息:\n");
            if IsBoundGlobal("AppendRetainedMirrorText") then
                ValueGlobal("AppendRetainedMirrorText")("        轨道信息:\n");
            fi;
            Print("          划分 (Partition): ", o.partition, "\n");
            if IsBoundGlobal("AppendRetainedMirrorText") then
                ValueGlobal("AppendRetainedMirrorText")(Concatenation("          划分 (Partition): ", String(o.partition), "\n"));
            fi;
            Print("          ab-diagram: ", o.ab_diagram, "\n");
            if IsBoundGlobal("AppendRetainedMirrorText") then
                ValueGlobal("AppendRetainedMirrorText")(Concatenation("          ab-diagram: ", String(o.ab_diagram), "\n"));
            fi;
            
            if subsystem_roots <> fail then
                # 将字符串列表转换为单行字符串（用换行符分隔）
                ab_str := "";
                for i in [1..Length(o.ab_diagram)] do
                    if i > 1 then ab_str := Concatenation(ab_str, "\n"); fi;
                    ab_str := Concatenation(ab_str, o.ab_diagram[i]);
                od;
                
                # 调用 Builder
                if IsBoundGlobal("GlobalAmbientType") and ValueGlobal("GlobalAmbientType") = "F4" then
                    Print("        [模块4重载] 重新加载 F4 专用 SL2 Builder\n");
                    if IsBoundGlobal("AppendRetainedMirrorText") then
                        ValueGlobal("AppendRetainedMirrorText")("        [模块4重载] 重新加载 F4 专用 SL2 Builder\n");
                    fi;
                    Read("../05_SL2_Triple_Construction/SL2_Triple_Builder.g");
                fi;
                if IsBoundGlobal("PrintSL2Triple") then
                    # 计算 ab-diagram 需要的边数 = sum(len(row)-1)
                    # 并将根列表扩展到足够长，避免出现占位符
                    total_edges := 0;
                    for row in o.ab_diagram do
                        if Length(row) > 1 then
                            total_edges := total_edges + (Length(row) - 1);
                        fi;
                    od;
                    n_base := Length(subsystem_roots);
                    if IsBoundGlobal("GlobalAmbientType") and ValueGlobal("GlobalAmbientType") = "F4" then
                        GlobalModule4CurrentRealFormInfo := rec(type := real_form_type, params := params);
                        ValueGlobal("PrintSL2Triple")(ab_str, real_form_type, params[1], subsystem_roots, fail);
                        if IsBoundGlobal("GlobalModule4CurrentRealFormInfo") then
                            Unbind(GlobalModule4CurrentRealFormInfo);
                        fi;
                    elif n_base > 0 and total_edges > n_base then
                        expanded_roots := [];
                        for i in [1..total_edges] do
                            Add(expanded_roots, subsystem_roots[((i - 1) mod n_base) + 1]);
                        od;
                        GlobalModule4CurrentRealFormInfo := rec(type := real_form_type, params := params);
                        ValueGlobal("PrintSL2Triple")(ab_str, real_form_type, params[1], expanded_roots, fail);
                        if IsBoundGlobal("GlobalModule4CurrentRealFormInfo") then
                            Unbind(GlobalModule4CurrentRealFormInfo);
                        fi;
                    else
                        GlobalModule4CurrentRealFormInfo := rec(type := real_form_type, params := params);
                        ValueGlobal("PrintSL2Triple")(ab_str, real_form_type, params[1], subsystem_roots, fail);
                        if IsBoundGlobal("GlobalModule4CurrentRealFormInfo") then
                            Unbind(GlobalModule4CurrentRealFormInfo);
                        fi;
                    fi;
                fi;
            fi;
            Print("\n");
            if IsBoundGlobal("AppendRetainedMirrorText") then
                ValueGlobal("AppendRetainedMirrorText")("\n");
            fi;
        od;
    fi;
end;;
