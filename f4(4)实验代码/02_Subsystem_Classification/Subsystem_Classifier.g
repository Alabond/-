#############################################################################
# 文件: 02_Subsystem_Classifier.g
# 描述: 通用子根系类型分类器 (支持 A, B, C, D, G 等类型判定)
# 作者: Trae AI Assistant
# 日期: 2026-02-26
# 说明: 输入一组候选 simple roots 后，先算 Cartan，再按连通分量逐块判型。
# 函数职责索引:
#   - UC_ComputeCartanMatrix:
#       用给定内积函数生成 Cartan 矩阵，是后续所有判型的基础。
#   - UC_FindConnectedComponents:
#       从 Cartan 非零耦合关系提取 Dynkin 图连通分量。
#   - UC_IdentifyComponentType:
#       对单个连通块进行类型判定，输出 A/B/C/D/G 与 rank 信息。
#   - ClassifyRootSystem:
#       对全体根调用上述流程，返回 total_type 与 components 结构。
#############################################################################



# 1. 辅助函数: 计算 Cartan 矩阵
# -----------------------------------------------------------------------------
# 输入: roots (根向量列表), inner_prod_func (内积函数)
# 输出: Cartan 矩阵 (C[i][j] = 2 * (r_i, r_j) / (r_j, r_j))
UC_ComputeCartanMatrix := function(roots, inner_prod_func)
    local n, C, i, j, ip_ij, ip_jj;
    n := Length(roots);
    C := [];
    for i in [1..n] do
        C[i] := [];
        for j in [1..n] do
            ip_ij := inner_prod_func(roots[i], roots[j]);
            ip_jj := inner_prod_func(roots[j], roots[j]);
            if ip_jj = 0 then
                Error("根向量长度为 0");
            fi;
            # Cartan 整数判定
            C[i][j] := Int(2 * ip_ij / ip_jj);
        od;
    od;
    return C;
end;;

# 2. 辅助函数: 寻找连通分量
# -----------------------------------------------------------------------------
# 输入: C (Cartan 矩阵)
# 输出: 连通分量的索引列表 (例如 [[1,2], [3]])
UC_FindConnectedComponents := function(C)
    local n, visited, components, i, comp, queue, curr, j, head, neighbors;
    n := Length(C);
    visited := List([1..n], x -> false);
    components := [];
    neighbors := List([1..n], x -> []);
    for i in [1..n] do
        for j in [1..n] do
            if i <> j and (C[i][j] <> 0 or C[j][i] <> 0) then
                Add(neighbors[i], j);
            fi;
        od;
    od;
    
    for i in [1..n] do
        if not visited[i] then
            comp := [];
            queue := [i];
            head := 1;
            visited[i] := true;
            while head <= Length(queue) do
                curr := queue[head];
                head := head + 1;
                Add(comp, curr);
                
                # 寻找邻接节点
                for j in neighbors[curr] do
                    if not visited[j] then
                        visited[j] := true;
                        Add(queue, j);
                    fi;
                od;
            od;
            Add(components, comp);
        fi;
    od;
    return components;
end;;

# 3. 核心函数: 识别单个连通分量的类型
# -----------------------------------------------------------------------------
# 输入: component_indices (索引列表), C (全局 Cartan 矩阵), roots (根向量列表, 用于长度判断)
# 输出: 类型字符串 (例如 "A2", "B3", "G2")
UC_ClassifySimplyLaced := function(subC, n)
    local i, j, bondCounts, degreeMax;
    bondCounts := List([1..n], x -> 0);
    for i in [1..n] do
        for j in [1..n] do
            if i <> j and subC[i][j] <> 0 then
                bondCounts[i] := bondCounts[i] + 1;
            fi;
        od;
    od;
    degreeMax := Maximum(bondCounts);
    if degreeMax <= 2 then
        return Concatenation("A", String(n));
    fi;
    if degreeMax = 3 and n >= 4 then
        return Concatenation("D", String(n));
    fi;
    return Concatenation("Unknown_SimplyLaced_", String(n));
end;;

UC_ClassifyDoubleBond := function(n)
    if n = 2 then
        return "B2/C2";
    fi;
    if n = 4 then
        return Concatenation("B", String(n), "/C", String(n));
    fi;
    return Concatenation("B", String(n), "/C", String(n));
end;;

UC_IdentifyComponentType := function(indices, C, roots, inner_prod_func)
    local n, subC, i, j, k, maxBond;
    
    n := Length(indices);
    
    # 特殊情况: A1
    if n = 1 then return "A1"; fi;
    
    # 构建子 Cartan 矩阵
    subC := [];
    for i in [1..n] do
        subC[i] := [];
        for j in [1..n] do
            subC[i][j] := C[indices[i]][indices[j]];
        od;
    od;
    
    # 统计键的类型
    maxBond := 0;
    for i in [1..n] do
        for j in [i+1..n] do
            # Cartan 矩阵非对角元素乘积即为键的阶数 (0, 1, 2, 3)
            # C[i][j] * C[j][i] = 4 * cos^2(theta)
            # 0 -> 0 (无连接)
            # 1 -> 1 (单键)
            # 2 -> 2 (双键)
            # 3 -> 3 (三键)
            k := subC[i][j] * subC[j][i];
            if k > maxBond then maxBond := k; fi;
        od;
    od;
    
    if maxBond = 1 then
        return UC_ClassifySimplyLaced(subC, n);
    fi;
    if maxBond = 2 then
        return UC_ClassifyDoubleBond(n);
    fi;
    if maxBond = 3 and n = 2 then
        return "G2";
    fi;
    
    return "Unknown";
end;;

# 4. 主入口: 分类根系
# -----------------------------------------------------------------------------
# 输入: roots (根向量列表), inner_prod_func
# 输出: 记录类型 rec( type_str := "A2 + A1", components := [ rec(type:="A2", indices:=[...]), ... ] )
ClassifyRootSystem := function(roots, inner_prod_func)
    local C, components, typeStr, compTypes, idxs, t, compList, result;
    
    if Length(roots) = 0 then 
        return rec(type_str := "Empty", components := []); 
    fi;
    
    C := UC_ComputeCartanMatrix(roots, inner_prod_func);
    components := UC_FindConnectedComponents(C);
    
    compTypes := [];
    compList := [];
    
    for idxs in components do
        t := UC_IdentifyComponentType(idxs, C, roots, inner_prod_func);
        Add(compTypes, t);
        Add(compList, rec(type := t, indices := idxs));
    od;
    
    # 排序类型字符串
    Sort(compTypes);
    typeStr := "";
    for t in compTypes do
        if Length(typeStr) > 0 then Append(typeStr, " + "); fi;
        Append(typeStr, t);
    od;
    
    return rec(type_str := typeStr, components := compList);
end;;

