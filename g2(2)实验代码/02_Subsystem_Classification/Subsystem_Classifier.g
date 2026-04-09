#############################################################################
# 文件: 02_Subsystem_Classifier.g
# 描述: 通用子根系类型分类器 (支持 A, B, C, D, G 等类型判定)
# 作者: Trae AI Assistant
# 日期: 2026-02-26
#############################################################################



# 1. 辅助函数: 计算 Cartan 矩阵
# -----------------------------------------------------------------------------
# 输入: roots (根向量列表), inner_prod_func (内积函数)
# 输出: Cartan 矩阵 (C[i][j] = 2 * (r_i, r_j) / (r_j, r_j))
# 按给定内积函数计算一组根的 Cartan 矩阵。
# 这是后续连通分量划分与 Dynkin 类型识别的基础数据。
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
# 根据 Cartan 矩阵中的非零耦合关系提取 Dynkin 图连通分量。
# 返回的是索引分组，而不是根本身，便于后续继续引用原始数据。
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
# 对 simply-laced 连通分量进行粗分类。
# 这里主要通过每个节点的度数模式区分 A_n 与 D_n。
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

# 对含双键的分量做类型标签归类。
# 当前实现只区分成 B_n/C_n 的并列表达，不进一步定向区分长短根方向。
UC_ClassifyDoubleBond := function(n)
    if n = 2 then
        return "B2/C2";
    fi;
    if n = 4 then
        return Concatenation("B", String(n), "/C", String(n));
    fi;
    return Concatenation("B", String(n), "/C", String(n));
end;;

# 识别单个连通分量的根系类型。
# 它先截取子 Cartan 矩阵，再根据最大键阶数与节点规模识别 A、D、B/C、G2 等类型。
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
# 对整组根做统一根系分类。
# 返回值同时包含整体类型字符串和每个连通分量的局部类型，供上游模块直接消费。
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
