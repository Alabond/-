#############################################################################
# 文件: SL2_Triple_Builder.g
# 描述: 从 ab-diagram 构造 sl2-triple (h, e, f)
# 参考: 项目内部实验流程与模块化约定
#############################################################################

# =============================================================================
# 基础辅助函数（先定义，避免未绑定）
# =============================================================================

# 浮点绝对值
# 统一计算数值的浮点绝对值。
# 该函数把有理数或整数先转换成浮点，再供误差阈值比较使用。
AbsFloat := function(x)
    local y;
    if IsFloat(x) then
        y := x;
    else
        y := Float(x);
    fi;
    if y >= 0.0 then return y; else return -y; fi;
end;;

# 记录 normal triple 校验统计信息。
# 每次模块 4 完成一次判定后都会更新总次数与通过次数。
RecordNormalTripleStats := function(is_pass)
    if not IsBoundGlobal("GlobalNormalTripleTotalChecks") then
        GlobalNormalTripleTotalChecks := 0;
    fi;
    if not IsBoundGlobal("GlobalNormalTriplePassChecks") then
        GlobalNormalTriplePassChecks := 0;
    fi;
    GlobalNormalTripleTotalChecks := GlobalNormalTripleTotalChecks + 1;
    if is_pass then
        GlobalNormalTriplePassChecks := GlobalNormalTriplePassChecks + 1;
    fi;
end;;

# =============================================================================
# 辅助函数
# =============================================================================

# 解析 ab-diagram 字符串为行列表
# 约定：
#   - 行以换行符分隔
#   - 每行长度 L 表示该 Jordan 块长度
# 示例：
#   "aba"     -> ["aba"]
#   "ab\na"   -> ["ab","a"]
# 将多行 ab-diagram 字符串拆分成行列表。
# 例如 "ab\nabab" 会被解析为 [ "ab", "abab" ]。
ParseABDiagram := function(diagram_str)
    local rows, row, ch;
    rows := [];
    row := "";
    for ch in diagram_str do
        if ch = '\n' then
            Add(rows, row);
            row := "";
        else
            row := Concatenation(row, [ch]);
        fi;
    od;
    if row <> "" then Add(rows, row); fi;
    return rows;
end;;

# 将 G2 根字符串解析为二维简单根坐标。
# 这里只处理 a1、a2 的线性组合，供 G2 专用链构造与校验使用。
ParseRootStringG2Simple := function(root_str)
    local v, coeff, idx, sign, num_str;
    v := [0, 0];
    root_str := Filtered(root_str, c -> c <> ' ');
    if root_str = "0" then return v; fi;
    idx := 1;
    while idx <= Length(root_str) do
        sign := 1;
        if root_str[idx] = '-' then
            sign := -1; idx := idx + 1;
        elif root_str[idx] = '+' then
            idx := idx + 1;
        fi;
        num_str := "";
        while idx <= Length(root_str) and root_str[idx] in "0123456789" do
            num_str := Concatenation(num_str, [root_str[idx]]);
            idx := idx + 1;
        od;
        if num_str = "" then coeff := 1; else coeff := Int(num_str); fi;
        coeff := coeff * sign;
        if idx <= Length(root_str) and root_str[idx] = 'a' then
            idx := idx + 1;
            if idx <= Length(root_str) and root_str[idx] in "12" then
                if root_str[idx] = '1' then v[1] := v[1] + coeff;
                elif root_str[idx] = '2' then v[2] := v[2] + coeff;
                fi;
                idx := idx + 1;
            fi;
        else
            break;
        fi;
    od;
    return v;
end;;

G2AmbientRoots := function()
    return [
        [1, 0], [0, 1], [1, 1], [2, 1], [3, 1], [3, 2],
        [-1, 0], [0, -1], [-1, -1], [-2, -1], [-3, -1], [-3, -2]
    ];
end;;

ConvertToG2RootVector := function(root_item)
    if IsString(root_item) then
        return ParseRootStringG2Simple(root_item);
    fi;
    if IsList(root_item) and Length(root_item) = 2 and IsInt(root_item[1]) and IsInt(root_item[2]) then
        return [root_item[1], root_item[2]];
    fi;
    return fail;
end;;

BuildGCJFromJ_G2 := function(subsystem_roots)
    local j_vecs, i, v, ambient, delta_vecs, det, n1, n2, key, keys, delta_names, h_basis, x_basis;
    j_vecs := [];
    if not IsList(subsystem_roots) then
        return fail;
    fi;
    for i in [1..Length(subsystem_roots)] do
        v := ConvertToG2RootVector(subsystem_roots[i]);
        if v <> fail then
            Add(j_vecs, v);
        fi;
    od;
    if Length(j_vecs) = 0 then
        return fail;
    fi;
    ambient := G2AmbientRoots();
    delta_vecs := [];
    keys := [];
    if Length(j_vecs) = 1 then
        for v in ambient do
            if j_vecs[1][1] * v[2] - j_vecs[1][2] * v[1] = 0 then
                key := FormatRoot(v);
                if Position(keys, key) = fail then
                    Add(keys, key);
                    Add(delta_vecs, v);
                fi;
            fi;
        od;
    else
        det := j_vecs[1][1] * j_vecs[2][2] - j_vecs[1][2] * j_vecs[2][1];
        if det = 0 then
            for v in ambient do
                if j_vecs[1][1] * v[2] - j_vecs[1][2] * v[1] = 0 then
                    key := FormatRoot(v);
                    if Position(keys, key) = fail then
                        Add(keys, key);
                        Add(delta_vecs, v);
                    fi;
                fi;
            od;
        else
            for v in ambient do
                n1 := v[1] * j_vecs[2][2] - v[2] * j_vecs[2][1];
                n2 := j_vecs[1][1] * v[2] - j_vecs[1][2] * v[1];
                if IsInt(n1 / det) and IsInt(n2 / det) then
                    key := FormatRoot(v);
                    if Position(keys, key) = fail then
                        Add(keys, key);
                        Add(delta_vecs, v);
                    fi;
                fi;
            od;
        fi;
    fi;
    delta_names := List(delta_vecs, r -> FormatRoot(r));
    Sort(delta_names);
    h_basis := [ "H_{a1}", "H_{a2}" ];
    x_basis := List(delta_names, r -> Concatenation("X_{[", r, "]}"));
    return rec(
        j_vectors := j_vecs,
        delta_vectors := delta_vecs,
        delta_roots := delta_names,
        cartan_basis := h_basis,
        x_basis := x_basis,
        dim_est := 2 + Length(delta_names)
    );
end;;

PrintGCJStep1 := function(gcj)
    if gcj = fail then
        Print("        [M4-STEP1] (g_c)_J 构造失败\n");
        return;
    fi;
    Print("        [M4-STEP1] (g_c)_J 构造\n");
    Print("        [M4-STEP1] J(向量) = ", gcj.j_vectors, "\n");
    Print("        [M4-STEP1] Delta_J = ", gcj.delta_roots, "\n");
    Print("        [M4-STEP1] h_c 基 = ", gcj.cartan_basis, "\n");
    Print("        [M4-STEP1] X_phi(phi∈Delta_J) = ", gcj.x_basis, "\n");
    Print("        [M4-STEP1] dim((g_c)_J) 估计 = ", gcj.dim_est, "\n");
end;;

ResolveModule4RealFormInfo := function(type, rank)
    local info;
    if IsBoundGlobal("GlobalModule4CurrentRealFormInfo") then
        info := ValueGlobal("GlobalModule4CurrentRealFormInfo");
        if IsRecord(info) and IsBound(info.type) and IsBound(info.params) then
            return info;
        fi;
    fi;
    if type = "Composite" and IsBoundGlobal("GlobalCompositeComponentData") and IsList(ValueGlobal("GlobalCompositeComponentData")) and Length(ValueGlobal("GlobalCompositeComponentData")) = 1 then
        info := ValueGlobal("GlobalCompositeComponentData")[1];
        if IsRecord(info) and IsBound(info.real_form_type) and IsBound(info.params) then
            return rec(type := info.real_form_type, params := info.params);
        fi;
    fi;
    if type = "su_pq" then
        if IsInt(rank) and rank > 0 then
            return rec(type := "su_pq", params := [rank, 1]);
        fi;
    fi;
    if type = "sl_R" then
        if IsInt(rank) and rank > 0 then
            return rec(type := "sl_R", params := [rank]);
        fi;
    fi;
    return fail;
end;;

BuildThetaAutomorphismFromRealForm := function(real_form_info)
    local rf_type, p, q, rank_a, auts;
    if real_form_info = fail or not IsRecord(real_form_info) then
        return fail;
    fi;
    if not (IsBoundGlobal("LoadPackage") and LoadPackage("sla") = true) then
        return fail;
    fi;
    if not IsBound(real_form_info.type) then
        return fail;
    fi;
    rf_type := real_form_info.type;
    if rf_type = "su_pq" and IsBound(real_form_info.params) and Length(real_form_info.params) >= 2 then
        p := real_form_info.params[1];
        q := real_form_info.params[2];
        if IsInt(p) and IsInt(q) and p > 0 and q > 0 then
            rank_a := p + q - 1;
            if rank_a >= 1 then
                auts := FiniteOrderInnerAutomorphisms("A", rank_a, 2);
                if IsList(auts) and Length(auts) > 0 then
                    return auts[1];
                fi;
            fi;
        fi;
    fi;
    if rf_type = "sl_R" and IsBound(real_form_info.params) and Length(real_form_info.params) >= 1 then
        rank_a := real_form_info.params[1];
        if IsInt(rank_a) and rank_a >= 1 then
            auts := FiniteOrderInnerAutomorphisms("A", rank_a, 2);
            if IsList(auts) and Length(auts) > 0 then
                return auts[1];
            fi;
        fi;
    fi;
    return fail;
end;;

ThetaCentralizerDimensionInL0 := function(theta, x)
    local g, l0_basis, L, BL, rows, b, br, coeffs, r;
    g := Grading(theta);
    if not IsList(g) or Length(g) = 0 then
        return fail;
    fi;
    l0_basis := g[1];
    if not IsList(l0_basis) then
        return fail;
    fi;
    if Length(l0_basis) = 0 then
        return 0;
    fi;
    L := Source(theta);
    BL := Basis(L);
    rows := [];
    for b in l0_basis do
        br := b * x;
        coeffs := Coefficients(BL, br);
        Add(rows, coeffs);
    od;
    if Length(rows) = 0 then
        return 0;
    fi;
    r := RankMat(rows);
    return Length(l0_basis) - r;
end;;

ThetaCenterIntersectionDimensionInL0 := function(theta)
    local g, l0_basis, L, lk, zlk;
    g := Grading(theta);
    if not IsList(g) or Length(g) = 0 then
        return fail;
    fi;
    l0_basis := g[1];
    if not IsList(l0_basis) then
        return fail;
    fi;
    L := Source(theta);
    lk := Subalgebra(L, l0_basis);
    zlk := Intersection(LieCenter(L), lk);
    return Dimension(zlk);
end;;

ComputeThetaOrbitTriplesForModule4 := function(real_form_info)
    local theta, orbs, triples, i, t, cdim, noticed, selected, kac, zdim, L, g, l0_basis, lk, zlk, e_subalg, c_full, c, zc, triple_subalg, c_triple_full, c_triple, triple_cdim, zc_dim, noticed_triples;
    if real_form_info = fail then
        return rec(ok := false, reason := "无实形式上下文");
    fi;
    theta := BuildThetaAutomorphismFromRealForm(real_form_info);
    if theta = fail then
        return rec(ok := false, reason := "无法构造theta");
    fi;
    if IsBoundGlobal("SetInfoLevel") and IsBoundGlobal("InfoSLA") then
        SetInfoLevel(InfoSLA, 1);
    fi;
    orbs := NilpotentOrbitsOfThetaRepresentation(theta);
    L := Source(theta);
    g := Grading(theta);
    if not IsList(g) or Length(g) = 0 or not IsList(g[1]) then
        return rec(ok := false, reason := "无法获取 l∩k 的分次信息");
    fi;
    l0_basis := g[1];
    lk := Subalgebra(L, l0_basis);
    zlk := Intersection(LieCenter(L), lk);
    zdim := Dimension(zlk);
    triples := [];
    for i in [1..Length(orbs)] do
        t := orbs[i];
        e_subalg := Subalgebra(L, [t[3]]);
        c_full := LieCentralizer(L, e_subalg);
        c := Intersection(c_full, lk);
        cdim := Dimension(c);
        zc := Intersection(c, zlk);
        zc_dim := Dimension(zc);
        triple_subalg := Subalgebra(L, [t[1], t[2], t[3]]);
        c_triple_full := LieCentralizer(L, triple_subalg);
        c_triple := Intersection(c_triple_full, lk);
        triple_cdim := Dimension(c_triple);
        noticed := triple_cdim = zdim;
        Add(triples, rec(
            index := i,
            y := t[1],
            h := t[2],
            x := t[3],
            centralizer_dim_lk := cdim,
            center_intersection_dim := zc_dim,
            triple_centralizer_dim_lk := triple_cdim,
            center_dim_lk := zdim,
            noticed := noticed
        ));
    od;
    selected := fail;
    noticed_triples := [];
    for t in triples do
        if t.noticed then
            Add(noticed_triples, t);
        fi;
    od;
    for t in triples do
        if t.noticed then
            selected := t;
            break;
        fi;
    od;
    if selected = fail and Length(triples) > 0 then
        selected := triples[1];
    fi;
    if HasKacDiagram(theta) then
        kac := KacDiagram(theta).weights;
    else
        kac := fail;
    fi;
    return rec(
        ok := true,
        real_form := real_form_info,
        theta := theta,
        kac_weights := kac,
        center_dim_lk := zdim,
        triples := triples,
        noticed_triples := noticed_triples,
        selected := selected
    );
end;;

MirrorRetainedBuilderText := function(text)
    if IsBoundGlobal("AppendRetainedMirrorText") then
        ValueGlobal("AppendRetainedMirrorText")(text);
    fi;
end;;
PrintThetaOrbitStep2 := function(theta_data, gcj)
    local t;
    if theta_data = fail or not IsRecord(theta_data) or not IsBound(theta_data.ok) or not theta_data.ok then
        if theta_data <> fail and IsRecord(theta_data) and IsBound(theta_data.reason) then
            Print("        [M4-STEP2] Theta-orbit 路径跳过: ", theta_data.reason, "\n");
            MirrorRetainedBuilderText(Concatenation("        [M4-STEP2] Theta-orbit 路径跳过: ", String(theta_data.reason), "\n"));
        else
            Print("        [M4-STEP2] Theta-orbit 路径跳过\n");
            MirrorRetainedBuilderText("        [M4-STEP2] Theta-orbit 路径跳过\n");
        fi;
        return;
    fi;
    Print("        [M4-STEP2] NilpotentOrbitsOfThetaRepresentation\n");
    MirrorRetainedBuilderText("        [M4-STEP2] NilpotentOrbitsOfThetaRepresentation\n");
    Print("        [M4-STEP2] RealForm = ", theta_data.real_form.type, ", Params = ", theta_data.real_form.params, "\n");
    MirrorRetainedBuilderText(Concatenation("        [M4-STEP2] RealForm = ", theta_data.real_form.type, ", Params = ", String(theta_data.real_form.params), "\n"));
    if IsBound(theta_data.kac_weights) and theta_data.kac_weights <> fail then
        Print("        [M4-STEP2] theta Kac weights = ", theta_data.kac_weights, "\n");
        MirrorRetainedBuilderText(Concatenation("        [M4-STEP2] theta Kac weights = ", String(theta_data.kac_weights), "\n"));
    fi;
    if IsBound(theta_data.center_dim_lk) and theta_data.center_dim_lk <> fail then
        Print("        [M4-STEP2] dim(z(l)∩k) = ", theta_data.center_dim_lk, "\n");
        MirrorRetainedBuilderText(Concatenation("        [M4-STEP2] dim(z(l)∩k) = ", String(theta_data.center_dim_lk), "\n"));
    fi;
    Print("        [M4-STEP2] 轨道数 = ", Length(theta_data.triples), "\n");
    MirrorRetainedBuilderText(Concatenation("        [M4-STEP2] 轨道数 = ", String(Length(theta_data.triples)), "\n"));
    for t in theta_data.triples do
        Print("        [M4-STEP2] orbit#", t.index, ": x=", t.x, ", h=", t.h, ", y=", t.y, "\n");
        MirrorRetainedBuilderText(Concatenation("        [M4-STEP2] orbit#", String(t.index), ": x=", String(t.x), ", h=", String(t.h), ", y=", String(t.y), "\n"));
        Print("        [M4-STEP2] orbit#", t.index, ": dim(C_e)= ", t.centralizer_dim_lk, ", dim(C_e∩Z)= ", t.center_intersection_dim, ", dim((l∩k)^[x,e,f])= ", t.triple_centralizer_dim_lk, ", dim(z(l)∩k)= ", t.center_dim_lk, ", noticed=", t.noticed, "\n");
        MirrorRetainedBuilderText(Concatenation("        [M4-STEP2] orbit#", String(t.index), ": dim(C_e)= ", String(t.centralizer_dim_lk), ", dim(C_e∩Z)= ", String(t.center_intersection_dim), ", dim((l∩k)^[x,e,f])= ", String(t.triple_centralizer_dim_lk), ", dim(z(l)∩k)= ", String(t.center_dim_lk), ", noticed=", String(t.noticed), "\n"));
    od;
    if theta_data.selected <> fail then
        Print("        [M4-STEP2] 选中 triple: orbit#", theta_data.selected.index, " (noticed=", theta_data.selected.noticed, ")\n");
        MirrorRetainedBuilderText(Concatenation("        [M4-STEP2] 选中 triple: orbit#", String(theta_data.selected.index), " (noticed=", String(theta_data.selected.noticed), ")\n"));
    fi;
    if gcj <> fail and IsRecord(gcj) and IsBound(gcj.dim_est) then
        Print("        [M4-STEP2] noticed 判定所用 l=(g_c)_J 维数估计 = ", gcj.dim_est, "\n");
        MirrorRetainedBuilderText(Concatenation("        [M4-STEP2] noticed 判定所用 l=(g_c)_J 维数估计 = ", String(gcj.dim_est), "\n"));
    fi;
end;;

NegateRootStringG2 := function(root_str)
    return FormatRoot(-ParseRootStringG2Simple(root_str));
end;;

AddRootStringsG2 := function(root_str1, root_str2)
    return FormatRoot(ParseRootStringG2Simple(root_str1) + ParseRootStringG2Simple(root_str2));
end;;

TranslateSelectedThetaTripleToRootTripleG2 := function(theta_data, ab_diagram_str, subsystem_roots, colors_opt)
    local selected, real_form, L, BL, x_coeffs, y_coeffs, slot_map, chain, e_map, f_map, h_map, slot, ex, fy, hh;
    if theta_data = fail or not IsRecord(theta_data) or not IsBound(theta_data.ok) or not theta_data.ok then
        return fail;
    fi;
    if not IsBound(theta_data.selected) or theta_data.selected = fail then
        return fail;
    fi;
    if not IsBound(theta_data.real_form) or theta_data.real_form = fail then
        return fail;
    fi;
    selected := theta_data.selected;
    real_form := theta_data.real_form;
    if not IsBound(theta_data.theta) or theta_data.theta = fail then
        return fail;
    fi;
    L := Source(theta_data.theta);
    BL := Basis(L);
    x_coeffs := Coefficients(BL, selected.x);
    y_coeffs := Coefficients(BL, selected.y);
    slot_map := fail;
    if real_form.type = "sl_R" and Length(subsystem_roots) >= 1 then
        slot_map := [
            rec(root := subsystem_roots[1], x_index := 2, y_index := 1),
            rec(root := NegateRootStringG2(subsystem_roots[1]), x_index := 1, y_index := 2)
        ];
    elif real_form.type = "su_pq" and Length(subsystem_roots) >= 2 and IsBound(real_form.params)
         and Length(real_form.params) >= 2 and real_form.params[1] = 2 and real_form.params[2] = 1 then
        chain := BuildAIIIChainForG2(ab_diagram_str, subsystem_roots, colors_opt);
        if chain = fail or not IsList(chain) or Length(chain) < 2 then
            return fail;
        fi;
        slot_map := [
            rec(root := chain[1], x_index := 3, y_index := 6),
            rec(root := chain[2], x_index := 4, y_index := 1),
            rec(root := NegateRootStringG2(subsystem_roots[2]), x_index := 2, y_index := 5)
        ];
    fi;
    if slot_map = fail then
        return fail;
    fi;
    e_map := [];
    f_map := [];
    h_map := [];
    for slot in slot_map do
        ex := x_coeffs[slot.x_index];
        fy := y_coeffs[slot.y_index];
        if ex <> 0 then
            Add(e_map, [slot.root, ex]);
        fi;
        if fy <> 0 then
            Add(f_map, [slot.root, fy]);
        fi;
        hh := ex * fy;
        if hh <> 0 then
            Add(h_map, [slot.root, hh]);
        fi;
    od;
    if Length(e_map) = 0 and Length(f_map) = 0 and Length(h_map) = 0 then
        return fail;
    fi;
    return rec(
        h := FormatExprFromMap(h_map, "H"),
        e := FormatExprFromMap(e_map, "E"),
        f := FormatExprFromMap(f_map, "F"),
        source := "theta_noticed"
    );
end;;

CollectThetaRootTripleCandidatesG2 := function(theta_data, ab_diagram_str, subsystem_roots, colors_opt)
    local candidates, t, td, one;
    candidates := [];
    if theta_data = fail or not IsRecord(theta_data) or not IsBound(theta_data.ok) or not theta_data.ok then
        return candidates;
    fi;
    if IsBound(theta_data.triples) and IsList(theta_data.triples) then
        for t in theta_data.triples do
            if IsBound(t.noticed) and t.noticed then
                td := ShallowCopy(theta_data);
                td.selected := t;
                one := TranslateSelectedThetaTripleToRootTripleG2(td, ab_diagram_str, subsystem_roots, colors_opt);
                if one <> fail then
                    Add(candidates, rec(triple := one, orbit_index := t.index));
                fi;
            fi;
        od;
    fi;
    if Length(candidates) = 0 and IsBound(theta_data.selected) and theta_data.selected <> fail then
        one := TranslateSelectedThetaTripleToRootTripleG2(theta_data, ab_diagram_str, subsystem_roots, colors_opt);
        if one <> fail then
            Add(candidates, rec(triple := one, orbit_index := theta_data.selected.index));
        fi;
    fi;
    return candidates;
end;;

# 检查 ab-diagram 是否满足项目中的 P' 条件。
# 该性质用于部分经典型线性构造分支的先验可行性判断。
HasPropertyPPrime := function(diagram_str)
    local rows, lens, i, j, li, lj, cnt;
    rows := ParseABDiagram(diagram_str);
    lens := List(rows, r -> Length(r));
    for i in [1..Length(lens)] do
        li := lens[i];
        if (li mod 2) = 0 then
            cnt := 0;
            for j in [1..Length(lens)] do
                lj := lens[j];
                if (lj mod 2) = 1 and lj < li then
                    cnt := cnt + 1;
                fi;
            od;
            if (cnt mod 2) <> 0 then
                return false;
            fi;
        fi;
    od;
    for i in [1..Length(lens)] do
        li := lens[i];
        if (li mod 2) = 1 then
            cnt := 0;
            for j in [1..Length(lens)] do
                lj := lens[j];
                if (lj mod 2) = 0 and lj > li then
                    cnt := cnt + 1;
                fi;
            od;
            if (cnt mod 2) <> 0 then
                return false;
            fi;
        fi;
    od;
    return true;
end;;

# 检查 ab-diagram 是否满足 P 条件。
# 主要用于 AIII 相关线性构造路径的必要条件检测。
HasPropertyP := function(diagram_str)
    local rows, lens, i, j, li, lj, cnt, k, la, lb;
    rows := ParseABDiagram(diagram_str);
    lens := List(rows, r -> Length(r));
    for i in [1..Length(lens)] do
        li := lens[i];
        if (li mod 2) = 0 then
            cnt := 0;
            for j in [1..Length(lens)] do
                lj := lens[j];
                if (lj mod 2) = 1 and lj < li then
                    cnt := cnt + 1;
                fi;
            od;
            if (cnt mod 2) <> 0 then
                return false;
            fi;
        fi;
    od;
    for i in [1..Length(lens)] do
        if (lens[i] mod 2) = 1 then
            for j in [1..Length(lens)] do
                if (lens[j] mod 2) = 1 and j <> i then
                    la := lens[i];
                    lb := lens[j];
                    cnt := 0;
                    for k in [1..Length(lens)] do
                        if (lens[k] mod 2) = 0 then
                            if (la < lb and la < lens[k] and lens[k] < lb) or (lb < la and lb < lens[k] and lens[k] < la) then
                                cnt := cnt + 1;
                            fi;
                        fi;
                    od;
                    if (cnt mod 2) <> 0 then
                        return false;
                    fi;
                fi;
            od;
        fi;
    od;
    return true;
end;;

# 统一入口：根据需要规范化根，再复用通用构造
# 规范化输入根列表，使其适合后续 normal triple 构造。
# 该函数会根据类型、秩和颜色信息调整顺序、过滤占位项并统一表示。
NormalizeRootsForNormalTriple := function(type, rank, roots, colors_opt)
    local i, has_compact, has_noncompact, compact_idx, non_compact_idx,
          vec_C, vec_NC, valid_g2_roots, is_root, r1_vec, r2_vec, r1_str, r2_str, sum_vec,
          M, len1, len2, ip, ambient, ParseRootSimple, ParseRoot, v1, v2, nc1, nc2;
    
    if not IsList(roots) or Length(roots) = 0 then
        return roots;
    fi;
    
    if rank = 2 and Length(roots) = 2 and IsList(colors_opt) and Length(colors_opt) = 2 then
        has_compact := false;
        has_noncompact := false;
        compact_idx := 0;
        non_compact_idx := 0;
        for i in [1..2] do
            if colors_opt[i] = 0 then
                has_compact := true; compact_idx := i;
            elif colors_opt[i] = 1 then
                has_noncompact := true; non_compact_idx := i;
            fi;
        od;
        if has_compact or has_noncompact then
            ambient := "G2";
            if IsBoundGlobal("GlobalAmbientType") then
                ambient := ValueGlobal("GlobalAmbientType");
            fi;
            ParseRootSimple := ParseRootStringG2Simple;
            ParseRoot := fail;
            if ambient = "G2" then
                ParseRoot := ParseRootSimple;
            elif IsBoundGlobal("ParseRootString") then
                ParseRoot := ValueGlobal("ParseRootString");
            else
                return roots;
            fi;
            if not (has_compact and has_noncompact) then
                if IsBoundGlobal("IsRootBlack") or IsBoundGlobal("IsCompact") then
                    has_compact := false;
                    has_noncompact := false;
                    compact_idx := 0;
                    non_compact_idx := 0;
                    for i in [1..2] do
                        if IsBoundGlobal("IsRootBlack") then
                            if ValueGlobal("IsRootBlack")(ParseRoot(roots[i])) then
                                has_noncompact := true; non_compact_idx := i;
                            else
                                has_compact := true; compact_idx := i;
                            fi;
                        elif IsBoundGlobal("IsCompact") then
                            if ValueGlobal("IsCompact")(ParseRoot(roots[i])) then
                                has_compact := true; compact_idx := i;
                            else
                                has_noncompact := true; non_compact_idx := i;
                            fi;
                        fi;
                    od;
                fi;
            fi;
            if not (has_compact and has_noncompact) then
                return roots;
            fi;
            if ambient = "G2" then
                if (IsBoundGlobal("IsRootBlack") or IsBoundGlobal("IsCompact")) then
                    v1 := ParseRoot(roots[1]);
                    v2 := ParseRoot(roots[2]);
                    if IsBoundGlobal("IsRootBlack") then
                        nc1 := ValueGlobal("IsRootBlack")(v1);
                        nc2 := ValueGlobal("IsRootBlack")(v2);
                    else
                        nc1 := not ValueGlobal("IsCompact")(v1);
                        nc2 := not ValueGlobal("IsCompact")(v2);
                    fi;
                    if nc1 <> nc2 then
                        if nc1 then
                            vec_NC := v1; vec_C := v2; non_compact_idx := 1; compact_idx := 2;
                        else
                            vec_NC := v2; vec_C := v1; non_compact_idx := 2; compact_idx := 1;
                        fi;
                        sum_vec := [vec_NC[1] + vec_C[1], vec_NC[2] + vec_C[2]];
                        valid_g2_roots := [
                            [1, 0], [0, 1], [1, 1], [2, 1], [3, 1], [3, 2],
                            [-1, 0], [0, -1], [-1, -1], [-2, -1], [-3, -1], [-3, -2]
                        ];
                        is_root := false;
                        for r1_vec in valid_g2_roots do
                            if r1_vec = sum_vec then
                                is_root := true;
                                break;
                            fi;
                        od;
                        if is_root then
                            if IsBoundGlobal("G2_InnerProdMatrix") then
                                M := ValueGlobal("G2_InnerProdMatrix");
                                ip := vec_NC[1] * M[1][1] * sum_vec[1] + vec_NC[1] * M[1][2] * sum_vec[2]
                                    + vec_NC[2] * M[2][1] * sum_vec[1] + vec_NC[2] * M[2][2] * sum_vec[2];
                                if ip > 0 then
                                    sum_vec := [ -sum_vec[1], -sum_vec[2] ];
                                fi;
                            fi;
                            r1_str := FormatRoot(vec_NC);
                            r2_str := FormatRoot(sum_vec);
                            return [ r1_str, r2_str ];
                        fi;
                    fi;
                fi;
            elif PositionSublist(type, "A") <> fail then
                return roots;
            fi;
            vec_NC := ParseRoot(roots[non_compact_idx]);
            vec_C := ParseRoot(roots[compact_idx]);
            if IsBoundGlobal("G2_InnerProdMatrix") then
                M := ValueGlobal("G2_InnerProdMatrix");
                len1 := vec_C[1] * M[1][1] * vec_C[1] + vec_C[1] * M[1][2] * vec_C[2] + vec_C[2] * M[2][1] * vec_C[1] + vec_C[2] * M[2][2] * vec_C[2];
                len2 := vec_NC[1] * M[1][1] * vec_NC[1] + vec_NC[1] * M[1][2] * vec_NC[2] + vec_NC[2] * M[2][1] * vec_NC[1] + vec_NC[2] * M[2][2] * vec_NC[2];
                ip := vec_C[1] * M[1][1] * vec_NC[1] + vec_C[1] * M[1][2] * vec_NC[2] + vec_C[2] * M[2][1] * vec_NC[1] + vec_C[2] * M[2][2] * vec_NC[2];
                if ip = 0 and ambient <> "G2" then
                    return roots;
                fi;
                if len1 = len2 and len1 <> 0 and ip * 2 = -len1 and ambient <> "G2" then
                    return roots;
                fi;
            fi;
            if IsList(vec_NC) and IsList(vec_C) and Length(vec_NC) = 2 and Length(vec_C) = 2 then
                sum_vec := [vec_NC[1] + vec_C[1], vec_NC[2] + vec_C[2]];
                valid_g2_roots := [
                    [1, 0], [0, 1], [1, 1], [2, 1], [3, 1], [3, 2],
                    [-1, 0], [0, -1], [-1, -1], [-2, -1], [-3, -1], [-3, -2]
                ];
                is_root := false;
                for r1_vec in valid_g2_roots do
                    if r1_vec = sum_vec then
                        is_root := true;
                        break;
                    fi;
                od;
                if is_root then
                    if IsBoundGlobal("G2_InnerProdMatrix") then
                        M := ValueGlobal("G2_InnerProdMatrix");
                        ip := vec_NC[1] * M[1][1] * sum_vec[1] + vec_NC[1] * M[1][2] * sum_vec[2]
                            + vec_NC[2] * M[2][1] * sum_vec[1] + vec_NC[2] * M[2][2] * sum_vec[2];
                        if ip > 0 then
                            sum_vec := [ -sum_vec[1], -sum_vec[2] ];
                        fi;
                    fi;
                    r1_str := FormatRoot(vec_NC);
                    r2_str := FormatRoot(sum_vec);
                    return [ r1_str, r2_str ];
                fi;
            fi;
        fi;
    fi;
    
    return roots;
end;;

# 将构造得到的 h 交给后续 K-orbit 输出模块。
# 当前 G2 路径主要通过它转发到 PrintKOrbitInfo。
DispatchKOrbitOutput := function(ambient, type, h_str, e_str, f_str)
    local effective_type;
    effective_type := type;
    if IsBoundGlobal("PrintKOrbitInfo") then
        if type = "Composite" and IsBoundGlobal("GlobalKOrbitSubsystemTypeOverride") then
            effective_type := ValueGlobal("GlobalKOrbitSubsystemTypeOverride");
        fi;
        GlobalCurrentSubsystemType := effective_type;
        ValueGlobal("PrintKOrbitInfo")(h_str);
    else
        PrintGenericHSummary(h_str);
    fi;
end;;

ApplySimpleReflectionG2Vec := function(v, alpha)
    local ip, coeff;
    ip := function(x, y)
        if IsBoundGlobal("InnerProduct") then
            return ValueGlobal("InnerProduct")(x, y);
        fi;
        return x[1] * 2 * y[1] + x[1] * (-3) * y[2] + x[2] * (-3) * y[1] + x[2] * 6 * y[2];
    end;
    coeff := 2 * ip(v, alpha) / ip(alpha, alpha);
    return [v[1] - coeff * alpha[1], v[2] - coeff * alpha[2]];
end;;

ApplyWKWordOnRootG2 := function(v, word, generators)
    local out, idx;
    out := [v[1], v[2]];
    for idx in word do
        out := ApplySimpleReflectionG2Vec(out, generators[idx]);
    od;
    return out;
end;;

SplitExprTermsByPlus := function(expr)
    local cleaned, terms, rest, pos;
    cleaned := Filtered(expr, c -> c <> '\n' and c <> '\r' and c <> '\t' and c <> '\\');
    while PositionSublist(cleaned, "  ") <> fail do
        cleaned := ReplacedString(cleaned, "  ", " ");
    od;
    if cleaned = "0" or cleaned = "" then
        return [];
    fi;
    terms := [];
    rest := cleaned;
    while true do
        pos := PositionSublist(rest, " + ");
        if pos = fail then
            Add(terms, rest);
            break;
        fi;
        Add(terms, rest{[1..pos-1]});
        if pos + 3 <= Length(rest) then
            rest := rest{[pos+3..Length(rest)]};
        else
            rest := "";
        fi;
    od;
    return terms;
end;;

BuildConjugacyKeyForExprG2 := function(expr, expected_kind, word, generators, parser)
    local terms, t, data, roots, v, wv;
    terms := SplitExprTermsByPlus(expr);
    roots := [];
    for t in terms do
        data := ParseTerm(t);
        if data.kind = expected_kind then
            v := parser(data.root);
            if v = fail then
                return fail;
            fi;
            wv := ApplyWKWordOnRootG2(v, word, generators);
            Add(roots, Concatenation(String(data.coeff), "@", FormatRoot(wv)));
        fi;
    od;
    return SortedList(roots);
end;;

IsWKConjugateFixingDeltaJG2 := function(t1, t2, subsystem_roots)
    local parser, gcj, delta_roots, delta_key, generators, words, stabilizer_words, word, moved, key1_h, key1_e, key1_f, key2_h, key2_e, key2_f;
    parser := ParseRootStringG2Simple;
    if IsBoundGlobal("ParseRootString") then
        parser := ValueGlobal("ParseRootString");
    fi;
    gcj := fail;
    if IsBoundGlobal("BuildGCJFromJ_G2") then
        gcj := ValueGlobal("BuildGCJFromJ_G2")(subsystem_roots);
    fi;
    if gcj = fail or not IsRecord(gcj) or not IsBound(gcj.delta_roots) then
        return false;
    fi;
    delta_roots := gcj.delta_roots;
    delta_key := SortedList(List(delta_roots, r -> FormatRoot(parser(r))));
    generators := [[1, 0], [-3, -2]];
    words := [[], [1], [2], [1, 2]];
    stabilizer_words := [];
    for word in words do
        moved := SortedList(List(delta_roots, r -> FormatRoot(ApplyWKWordOnRootG2(parser(r), word, generators))));
        if moved = delta_key then
            Add(stabilizer_words, word);
        fi;
    od;
    if Length(stabilizer_words) = 0 then
        return false;
    fi;
    key2_h := BuildConjugacyKeyForExprG2(t2.h, "H", [], generators, parser);
    key2_e := BuildConjugacyKeyForExprG2(t2.e, "E", [], generators, parser);
    key2_f := BuildConjugacyKeyForExprG2(t2.f, "F", [], generators, parser);
    if key2_h = fail or key2_e = fail or key2_f = fail then
        return false;
    fi;
    for word in stabilizer_words do
        key1_h := BuildConjugacyKeyForExprG2(t1.h, "H", word, generators, parser);
        key1_e := BuildConjugacyKeyForExprG2(t1.e, "E", word, generators, parser);
        key1_f := BuildConjugacyKeyForExprG2(t1.f, "F", word, generators, parser);
        if key1_h <> fail and key1_e <> fail and key1_f <> fail and key1_h = key2_h and key1_e = key2_e and key1_f = key2_f then
            return true;
        fi;
    od;
    return false;
end;;

FilterCompositeTriplesByWKConjugacyFixingDeltaJG2 := function(triples, subsystem_roots)
    local kept, i, j, is_dup;
    kept := [];
    for i in [1..Length(triples)] do
        is_dup := false;
        for j in [1..Length(kept)] do
            if IsWKConjugateFixingDeltaJG2(triples[i], kept[j], subsystem_roots) then
                is_dup := true;
                break;
            fi;
        od;
        if not is_dup then
            Add(kept, triples[i]);
        fi;
    od;
    return kept;
end;;

# 模块 4 的统一主入口。
# 负责根列表规范化、线性构造/组件构造、EFH 与 P/K 校验以及 K-orbit 联动输出。
PrintSL2TripleUnified := function(ab_diagram_str, type, rank, subsystem_roots, colors_opt)
    local normalized_roots, e_str, h_str, f_str, effective_colors, triple, efh_ok, pk_ok, diag_ok, gcj_data, real_form_info, theta_data, is_composite, theta_triple, k, m5_item, m5_count, raw_m5_count;

    normalized_roots := NormalizeRootsForNormalTriple(type, rank, subsystem_roots, colors_opt);
    effective_colors := colors_opt;
    is_composite := (type = "Composite" and IsBoundGlobal("GlobalCompositeComponentData") and IsList(ValueGlobal("GlobalCompositeComponentData")) and Length(ValueGlobal("GlobalCompositeComponentData")) > 0);
    if is_composite then
        gcj_data := fail;
        theta_data := fail;
        GlobalCurrentGCJData := fail;
        GlobalCurrentThetaOrbitData := fail;
        triple := BuildCompositeNormalTriple_FromComponents(ValueGlobal("GlobalCompositeComponentData"));
        Print("        [COMPOSITE-RAW] h = ", triple.h, "\n");
        MirrorRetainedBuilderText(Concatenation("        [COMPOSITE-RAW] h = ", String(triple.h), "\n"));
        Print("        [COMPOSITE-RAW] e = ", triple.e, "\n");
        MirrorRetainedBuilderText(Concatenation("        [COMPOSITE-RAW] e = ", String(triple.e), "\n"));
        Print("        [COMPOSITE-RAW] f = ", triple.f, "\n");
        MirrorRetainedBuilderText(Concatenation("        [COMPOSITE-RAW] f = ", String(triple.f), "\n"));
        Print("        [COMPOSITE] 按组件构造并求和完成\n");
        MirrorRetainedBuilderText("        [COMPOSITE] 按组件构造并求和完成\n");
    else
        gcj_data := BuildGCJFromJ_G2(subsystem_roots);
        GlobalCurrentGCJData := gcj_data;
        PrintGCJStep1(gcj_data);
        real_form_info := ResolveModule4RealFormInfo(type, rank);
        theta_data := ComputeThetaOrbitTriplesForModule4(real_form_info);
        GlobalCurrentThetaOrbitData := theta_data;
        PrintThetaOrbitStep2(theta_data, gcj_data);
        theta_triple := TranslateSelectedThetaTripleToRootTripleG2(theta_data, ab_diagram_str, normalized_roots, effective_colors);
        if theta_triple <> fail then
            triple := theta_triple;
            Print("        [M4-STEP3] noticed triple 已翻译到根形式\n");
            MirrorRetainedBuilderText("        [M4-STEP3] noticed triple 已翻译到根形式\n");
            Print("        [THETA-CANDIDATE] h = ", triple.h, "\n");
            MirrorRetainedBuilderText(Concatenation("        [THETA-CANDIDATE] h = ", String(triple.h), "\n"));
            Print("        [THETA-CANDIDATE] e = ", triple.e, "\n");
            MirrorRetainedBuilderText(Concatenation("        [THETA-CANDIDATE] e = ", String(triple.e), "\n"));
            Print("        [THETA-CANDIDATE] f = ", triple.f, "\n");
            MirrorRetainedBuilderText(Concatenation("        [THETA-CANDIDATE] f = ", String(triple.f), "\n"));
        else
            triple := BuildBDINormalTriple_Linear(ab_diagram_str, type, rank, normalized_roots, effective_colors);
        fi;
    fi;
    h_str := triple.h;
    e_str := triple.e;
    f_str := triple.f;
    Print("        >>> 模块 4 输出: SL2-Triple 构造 (Candidate)\n");
    MirrorRetainedBuilderText("        >>> 模块 4 输出: SL2-Triple 构造 (Candidate)\n");
    Print("        SL2-Triple (h, e, f):\n");
    MirrorRetainedBuilderText("        SL2-Triple (h, e, f):\n");
    Print("          h = ", h_str, "\n");
    MirrorRetainedBuilderText(Concatenation("          h = ", String(h_str), "\n"));
    Print("          e = ", e_str, "\n");
    MirrorRetainedBuilderText(Concatenation("          e = ", String(e_str), "\n"));
    Print("          f = ", f_str, "\n");
    MirrorRetainedBuilderText(Concatenation("          f = ", String(f_str), "\n"));

    Print("        [强制校验] 开始 VerifyEFH_NormalConditions…\n");
    MirrorRetainedBuilderText("        [强制校验] 开始 VerifyEFH_NormalConditions…\n");
    if IsBoundGlobal("VerifyEFH_NormalConditions") then
        efh_ok := ValueGlobal("VerifyEFH_NormalConditions")(h_str, e_str, f_str);
    else
        efh_ok := false;
        Print("        [强制校验] VerifyEFH_NormalConditions 未加载！\n");
        MirrorRetainedBuilderText("        [强制校验] VerifyEFH_NormalConditions 未加载！\n");
    fi;
    pk_ok := true;
    if IsBoundGlobal("VerifyPK_Subsystem_G2") then
        pk_ok := ValueGlobal("VerifyPK_Subsystem_G2")(h_str, e_str, f_str, normalized_roots, effective_colors);
        if pk_ok then
            Print("        [诊断] P/K 检查: PASS\n");
            MirrorRetainedBuilderText("        [诊断] P/K 检查: PASS\n");
        else
            Print("        [诊断] P/K 检查: FAIL\n");
            MirrorRetainedBuilderText("        [诊断] P/K 检查: FAIL\n");
        fi;
    fi;
    diag_ok := (efh_ok and pk_ok);
    if diag_ok then
        Print("        模块4判定: Normal\n");
        MirrorRetainedBuilderText("        模块4判定: Normal\n");
        Print("NORMAL_TRIPLE_CHECK: PASS\n");
        MirrorRetainedBuilderText("NORMAL_TRIPLE_CHECK: PASS\n");
        RecordNormalTripleStats(true);
    else
        Print("        模块4判定: Non-normal\n");
        MirrorRetainedBuilderText("        模块4判定: Non-normal\n");
        Print("NORMAL_TRIPLE_CHECK: FAIL\n");
        MirrorRetainedBuilderText("NORMAL_TRIPLE_CHECK: FAIL\n");
        RecordNormalTripleStats(false);
    fi;

    Print("\n");
    MirrorRetainedBuilderText("\n");
    if not IsBoundGlobal("PrintKOrbitInfo") and IsBoundGlobal("Read") then
        Read("../04_Orbit_Analysis/K_Orbit_Classifier.g");
    fi;
    if is_composite and IsRecord(triple) and IsBound(triple.all_triples) and IsList(triple.all_triples) and Length(triple.all_triples) > 0 then
        raw_m5_count := Length(triple.all_triples);
        triple.all_triples := FilterCompositeTriplesByWKConjugacyFixingDeltaJG2(triple.all_triples, subsystem_roots);
        if Length(triple.all_triples) < raw_m5_count then
            Print("        [M4-WK-FILTER] 按“W_k 共轭且固定 Delta_J”过滤: ", raw_m5_count, " -> ", Length(triple.all_triples), "\n");
            MirrorRetainedBuilderText(Concatenation("        [M4-WK-FILTER] 按“W_k 共轭且固定 Delta_J”过滤: ", String(raw_m5_count), " -> ", String(Length(triple.all_triples)), "\n"));
        else
            Print("        [M4-WK-FILTER] 无可过滤的 W_k 共轭 triple (固定 Delta_J)\n");
            MirrorRetainedBuilderText("        [M4-WK-FILTER] 无可过滤的 W_k 共轭 triple (固定 Delta_J)\n");
        fi;
        m5_count := Length(triple.all_triples);
        Print("        [M5-PROPAGATE] 向模块5传递组合 triple 数 = ", m5_count, "\n");
        MirrorRetainedBuilderText(Concatenation("        [M5-PROPAGATE] 向模块5传递组合 triple 数 = ", String(m5_count), "\n"));
        for k in [1..m5_count] do
            m5_item := triple.all_triples[k];
            Print("        [M5-PROPAGATE#", k, "] h = ", m5_item.h, "\n");
            MirrorRetainedBuilderText(Concatenation("        [M5-PROPAGATE#", String(k), "] h = ", String(m5_item.h), "\n"));
            Print("        [M5-PROPAGATE#", k, "] e = ", m5_item.e, "\n");
            MirrorRetainedBuilderText(Concatenation("        [M5-PROPAGATE#", String(k), "] e = ", String(m5_item.e), "\n"));
            Print("        [M5-PROPAGATE#", k, "] f = ", m5_item.f, "\n");
            MirrorRetainedBuilderText(Concatenation("        [M5-PROPAGATE#", String(k), "] f = ", String(m5_item.f), "\n"));
            if IsBound(m5_item.component_orbit_choice) then
                Print("        [M5-PROPAGATE#", k, "] 组件orbit选择 = ", m5_item.component_orbit_choice, "\n");
                MirrorRetainedBuilderText(Concatenation("        [M5-PROPAGATE#", String(k), "] 组件orbit选择 = ", String(m5_item.component_orbit_choice), "\n"));
            fi;
            DispatchKOrbitOutput("G2", type, m5_item.h, m5_item.e, m5_item.f);
        od;
    else
        DispatchKOrbitOutput("G2", type, h_str, e_str, f_str);
    fi;
end;;

# =============================================================================
# 核心构造逻辑
# =============================================================================

# ConstructEElement: 构造幂零元 e
# 输入: ab-diagram (字符串列表), 子系统类型 (如 "A2")
# 输出: 字符串描述，形式如 "E_{alpha1} + E_{alpha2} + ..."

# 映射通用根名称到 G2 子系统中的具体根
# 输入: 通用索引 (1, 2, ...), G2子系统根列表 (w_alpha2, neg_w_theta)
#
# 对于 A2 子系统 (秩 2)，有两个简单根:
# root 1: w_alpha2 (长根)
# root 2: neg_w_theta (扩展根，短根? 不，neg_w_theta 是 -3a1-2a2，长根?
# 在 G2 中:
# alpha1 (短), alpha2 (长)
# theta = 3a1 + 2a2 (长)
# 所以 w(theta) 是长根。
# 
# 我们的子系统是由 { w(alpha2), -w(theta) } 生成的。
# 这两个根都是长根。
# 它们构成的确实是 A2 (因为角度 120度，长度相等)。
#
# 因此，当我们在通用逻辑中说 "a1", "a2" 时，
# 我们需要将它们映射回 G2 的这两个具体根。
# 映射顺序:
# 通常我们按某种顺序排列简单根。
# 这里简单起见，我们定义:
#   Generic a1 -> Subsystem Root 1 (w_alpha2)
#   Generic a2 -> Subsystem Root 2 (neg_w_theta)
# (或者反过来，对于 A2 是对称的)



# ConstructHElement: 构造半单元 h
# 输入: ab-diagram (字符串), 子系统根列表 (字符串形式)
# 输出: "H_{...} + ..." 形式的字符串
# 原理: 单个长度为 L 的行，对应特征值序列 L-1, L-3, ..., -(L-1)。
# 将这些特征值在行内做前缀和，得到各简单根对应的系数 c_k，
# =============================================================================
# 主接口：打印 (h, e, f) 并联动 K-Orbit 分类输出
# =============================================================================

# 在未加载 K-orbit 分类器时输出 h 的兜底摘要。
# 它只做轻量解析，不参与真正的轨道分类。
PrintGenericHSummary := function(h_expr)
    local parts, i, t, coeffs, c, ch, num, sign, k, term, coef, inside;
    if IsBoundGlobal("PrintKOrbitInfo") then
        GlobalCurrentSubsystemType := "Unknown";
        ValueGlobal("PrintKOrbitInfo")(h_expr);
        return;
    fi;
    parts := [];
    term := "";
    for ch in h_expr do
        if ch = '+' then
            if term <> "" then Add(parts, term); fi;
            term := "";
        else
            term := Concatenation(term, [ch]);
        fi;
    od;
    if term <> "" then Add(parts, term); fi;
    coeffs := [];
    for t in parts do
        i := 1;
        while i <= Length(t) and (t[i] = ' ' or t[i] = '\t') do i := i + 1; od;
        sign := 1;
        if i <= Length(t) and t[i] = '-' then
            sign := -1; i := i + 1;
        fi;
        num := "";
        while i <= Length(t) and t[i] in "0123456789" do
            num := Concatenation(num, [t[i]]);
            i := i + 1;
        od;
        if num = "" then
            if i <= Length(t) and t[i] = 'H' then
                coef := 1 * sign;
            else
                coef := 0;
            fi;
        else
            coef := Int(num) * sign;
        fi;
        Add(coeffs, coef);
    od;
    Print("        >>> 模块 5 输出: H-摘要\n");
    Print("        H 项数: ", Length(coeffs), "\n");
end;;

# 合并两个“根名 → 系数”映射列表。
# 若同一根在两张映射表中都出现，则自动相加系数。
MergeTermMaps := function(map1, map2)
    local result, pair, found, p;
    result := ShallowCopy(map1);
    for pair in map2 do
        found := false;
        for p in result do
            if p[1] = pair[1] then
                p[2] := p[2] + pair[2];
                found := true;
                break;
            fi;
        od;
        if not found then
            Add(result, [pair[1], pair[2]]);
        fi;
    od;
    return result;
end;;

# 把系数映射重新格式化为 H/E/F 表达式字符串。
# 该函数是组件求和后重建表达式的统一出口。
FormatExprFromMap := function(map, kind)
    local expr, pair, coeff, root;
    expr := "";
    for pair in map do
        coeff := pair[2];
        root := pair[1];
        if expr <> "" then expr := Concatenation(expr, " + "); fi;
        if kind = "H" then
            if coeff = 1 then
                expr := Concatenation(expr, "H_{", root, "}");
            else
                expr := Concatenation(expr, String(coeff), "H_{", root, "}");
            fi;
        elif kind = "E" then
            if coeff = 1 then
                expr := Concatenation(expr, "E_{[", root, "]}");
            else
                expr := Concatenation(expr, String(coeff), "E_{[", root, "]}");
            fi;
        else
            if coeff = 1 then
                expr := Concatenation(expr, "F_{[", root, "]}");
            else
                expr := Concatenation(expr, String(coeff), "F_{[", root, "]}");
            fi;
        fi;
    od;
    if expr = "" then
        expr := "0";
    fi;
    return expr;
end;;

# 从多个已拆分的组件数据分别构造三元组并求和。
# 主要用于复合子系统，把每个组件的局部 normal triple 汇总成全局结果。
BuildCompositeNormalTriple_FromComponents := function(component_data)
    local comp, rows, i, comp_ab_str, comp_roots, comp_colors, comp_type, comp_rank, comp_triple,
          h_expr, e_expr, f_expr, normalized_comp_roots, comp_efh_ok, comp_pk_ok, comp_diag_ok,
          gcj_data, comp_real_form_info, theta_data, comp_idx, comp_candidates, candidate, cand_idx,
          all_component_candidates, combined, new_combined, base, item, merged, path, k;
    h_expr := "";
    e_expr := "";
    f_expr := "";
    all_component_candidates := [];
    comp_idx := 0;
    for comp in component_data do
        comp_idx := comp_idx + 1;
        rows := comp.ab_diagram;
        if IsString(rows) then
            comp_ab_str := rows;
        else
            comp_ab_str := "";
            for i in [1..Length(rows)] do
                if i > 1 then
                    comp_ab_str := Concatenation(comp_ab_str, "\n");
                fi;
                comp_ab_str := Concatenation(comp_ab_str, rows[i]);
            od;
        fi;
        comp_roots := comp.roots;
        comp_colors := fail;
        if IsBound(comp.colors) then
            comp_colors := comp.colors;
        fi;
        comp_type := "CompositeComponent";
        if IsBound(comp.real_form_type) then
            comp_type := comp.real_form_type;
        fi;
        comp_rank := Length(comp_roots);
        Print("        [COMP#", comp_idx, "] 开始模块4组件流程\n");
        gcj_data := BuildGCJFromJ_G2(comp_roots);
        GlobalCurrentGCJData := gcj_data;
        PrintGCJStep1(gcj_data);
        comp_real_form_info := fail;
        if IsBound(comp.real_form_type) and IsBound(comp.params) then
            comp_real_form_info := rec(type := comp.real_form_type, params := comp.params);
        fi;
        theta_data := ComputeThetaOrbitTriplesForModule4(comp_real_form_info);
        GlobalCurrentThetaOrbitData := theta_data;
        PrintThetaOrbitStep2(theta_data, gcj_data);
        normalized_comp_roots := NormalizeRootsForNormalTriple(comp_type, comp_rank, comp_roots, comp_colors);
        if IsBoundGlobal("GlobalAmbientType") and ValueGlobal("GlobalAmbientType") = "G2" and (comp_type = "su_pq" or PositionSublist(comp_type, "su(") <> fail or PositionSublist(comp_type, "AIII") <> fail) then
            normalized_comp_roots := comp_roots;
        fi;
        comp_candidates := CollectThetaRootTripleCandidatesG2(theta_data, comp_ab_str, normalized_comp_roots, comp_colors);
        if Length(comp_candidates) = 0 then
            if comp_rank = 1 and Length(normalized_comp_roots) >= 1 then
                comp_triple := rec(
                    h := Concatenation("H_{", normalized_comp_roots[1], "}"),
                    e := Concatenation("E_{[", normalized_comp_roots[1], "]}"),
                    f := Concatenation("F_{[", normalized_comp_roots[1], "]}")
                );
            else
                comp_triple := BuildBDINormalTriple_Linear(comp_ab_str, comp_type, comp_rank, normalized_comp_roots, comp_colors);
            fi;
            comp_candidates := [rec(triple := comp_triple, orbit_index := fail)];
        else
            Print("        [组件求和输入] 来源: noticed triple 根形式翻译\n");
            Print("        [COMP#", comp_idx, "] noticed 候选数 = ", Length(comp_candidates), "\n");
        fi;
        cand_idx := 0;
        for candidate in comp_candidates do
            cand_idx := cand_idx + 1;
            Print("        [组件候选#", cand_idx, "] h = ", candidate.triple.h, "\n");
            Print("        [组件候选#", cand_idx, "] e = ", candidate.triple.e, "\n");
            Print("        [组件候选#", cand_idx, "] f = ", candidate.triple.f, "\n");
            if candidate.orbit_index <> fail then
                Print("        [组件候选#", cand_idx, "] 来自 orbit#", candidate.orbit_index, "\n");
            fi;
        od;
        comp_triple := comp_candidates[1].triple;
        Print("        [组件诊断] 类型=", comp_type, ", 根=", comp_roots, "\n");
        Print("        [组件诊断] 开始 VerifyEFH_NormalConditions…\n");
        comp_efh_ok := false;
        if IsBoundGlobal("VerifyEFH_NormalConditions") then
            comp_efh_ok := ValueGlobal("VerifyEFH_NormalConditions")(comp_triple.h, comp_triple.e, comp_triple.f);
        fi;
        comp_pk_ok := true;
        if IsBoundGlobal("VerifyPK_Subsystem_G2") then
            comp_pk_ok := ValueGlobal("VerifyPK_Subsystem_G2")(comp_triple.h, comp_triple.e, comp_triple.f, normalized_comp_roots, comp_colors);
            if comp_pk_ok then
                Print("        [组件诊断] P/K 检查: PASS\n");
            else
                Print("        [组件诊断] P/K 检查: FAIL\n");
            fi;
        fi;
        comp_diag_ok := (comp_efh_ok and comp_pk_ok);
        if comp_diag_ok then
            Print("        [组件诊断] NORMAL_TRIPLE_CHECK: PASS\n");
        else
            Print("        [组件诊断] NORMAL_TRIPLE_CHECK: FAIL\n");
        fi;
        Add(all_component_candidates, comp_candidates);
    od;
    combined := [rec(h := "", e := "", f := "", path := [])];
    for comp_candidates in all_component_candidates do
        new_combined := [];
        for base in combined do
            for item in comp_candidates do
                merged := rec(h := base.h, e := base.e, f := base.f);
                if item.triple.h <> "0" then
                    if merged.h <> "" then merged.h := Concatenation(merged.h, " + "); fi;
                    merged.h := Concatenation(merged.h, item.triple.h);
                fi;
                if item.triple.e <> "0" then
                    if merged.e <> "" then merged.e := Concatenation(merged.e, " + "); fi;
                    merged.e := Concatenation(merged.e, item.triple.e);
                fi;
                if item.triple.f <> "0" then
                    if merged.f <> "" then merged.f := Concatenation(merged.f, " + "); fi;
                    merged.f := Concatenation(merged.f, item.triple.f);
                fi;
                path := ShallowCopy(base.path);
                Add(path, item.orbit_index);
                Add(new_combined, rec(h := merged.h, e := merged.e, f := merged.f, path := path));
            od;
        od;
        combined := new_combined;
    od;
    Print("        [COMPOSITE] 组件候选笛卡尔组合数 = ", Length(combined), "\n");
    for k in [1..Length(combined)] do
        if combined[k].h = "" then combined[k].h := "0"; fi;
        if combined[k].e = "" then combined[k].e := "0"; fi;
        if combined[k].f = "" then combined[k].f := "0"; fi;
        Print("        [COMPOSITE-RAW#", k, "] h = ", combined[k].h, "\n");
        Print("        [COMPOSITE-RAW#", k, "] e = ", combined[k].e, "\n");
        Print("        [COMPOSITE-RAW#", k, "] f = ", combined[k].f, "\n");
        Print("        [COMPOSITE-RAW#", k, "] 组件orbit选择 = ", combined[k].path, "\n");
    od;
    h_expr := combined[1].h;
    e_expr := combined[1].e;
    f_expr := combined[1].f;
    return rec(
        h := h_expr,
        e := e_expr,
        f := f_expr,
        all_triples := List(combined, c -> rec(h := c.h, e := c.e, f := c.f, component_orbit_choice := c.path))
    );
end;;

# 把 sl2-triple 的系数归一到单位尺度。
# 用于后续重整、比较或校验时消除整体倍数差异。
NormalizeSL2Triple_Unit := function(h_str, e_str, f_str)
    local e_terms, f_terms, i, parts, pr, pair, roots, seen, r, e_norm, f_norm, h_norm;
    e_terms := [];
    f_terms := [];
    if e_str <> "0" then
        for pr in SplitString(e_str, " + ") do
            pair := ValueGlobal("ParseTerm")(pr);
            if pair.kind = "E" then Add(e_terms, [pair.root, pair.coeff]); fi;
        od;
    fi;
    if f_str <> "0" then
        for pr in SplitString(f_str, " + ") do
            pair := ValueGlobal("ParseTerm")(pr);
            if pair.kind = "F" then Add(f_terms, [pair.root, pair.coeff]); fi;
        od;
    fi;
    if Length(e_terms) <> Length(f_terms) then
        return rec(ok := false, h := h_str, e := e_str, f := f_str);
    fi;
    for i in [1..Length(e_terms)] do
        if e_terms[i][1] <> f_terms[i][1] then
            return rec(ok := false, h := h_str, e := e_str, f := f_str);
        fi;
    od;
    roots := [];
    seen := [];
    for i in [1..Length(e_terms)] do
        r := e_terms[i][1];
        if not r in seen then
            Add(seen, r);
            Add(roots, r);
        fi;
    od;
    e_norm := "";
    f_norm := "";
    h_norm := "";
    for i in [1..Length(roots)] do
        r := roots[i];
        if e_norm <> "" then e_norm := Concatenation(e_norm, " + "); fi;
        e_norm := Concatenation(e_norm, "E_{[", r, "]}");
        if f_norm <> "" then f_norm := Concatenation(f_norm, " + "); fi;
        f_norm := Concatenation(f_norm, "F_{[", r, "]}");
        if h_norm <> "" then h_norm := Concatenation(h_norm, " + "); fi;
        h_norm := Concatenation(h_norm, "H_{", r, "}");
    od;
    if e_norm = "" then e_norm := "0"; fi;
    if f_norm = "" then f_norm := "0"; fi;
    if h_norm = "" then h_norm := "0"; fi;
    return rec(ok := true, h := h_norm, e := e_norm, f := f_norm);
end;;

# 求解有理系数线性方程组。
# 用于经典型线性构造中精确求系数，尽量避免浮点误差。
SolveRationalLinearSystem := function(A, b)
    local n, Ab, i, j, k, pivot, factor, temp, rank, x;
    n := Length(A);
    if n = 0 or Length(b) <> n then return fail; fi;
    # 构造增广矩阵 [A|b]
    Ab := [];
    for i in [1..n] do
        Ab[i] := ShallowCopy(A[i]);
        Add(Ab[i], b[i]);
    od;
    # 前向消元
    rank := 0;
    for k in [1..n] do
        # 找主元
        pivot := 0;
        for i in [k..n] do
            if Ab[i][k] <> 0 then pivot := i; break; fi;
        od;
        if pivot = 0 then continue; fi;
        rank := rank + 1;
        # 交换行
        if pivot <> k then
            temp := Ab[k];
            Ab[k] := Ab[pivot];
            Ab[pivot] := temp;
        fi;
        # 消元
        for i in [k+1..n] do
            factor := Ab[i][k] / Ab[k][k];
            for j in [k..n+1] do
                Ab[i][j] := Ab[i][j] - factor * Ab[k][j];
            od;
        od;
    od;
    # 回代
    x := List([1..n], i -> 0);
    for i in [n, n-1..1] do
        if Ab[i][i] = 0 then
            if Ab[i][n+1] <> 0 then return fail; fi;
            x[i] := 0;  # 自由变量，这里简单设为 0
        else
            x[i] := Ab[i][n+1];
            for j in [i+1..n] do
                x[i] := x[i] - Ab[i][j] * x[j];
            od;
            x[i] := x[i] / Ab[i][i];
        fi;
    od;
    return x;
end;;

# 判断一个二维向量是否是当前实现支持的 G2 根。
# 主要用于 G2 特化链构造与输入防御检查。
IsKnownG2RootVec := function(v)
    local roots, r;
    roots := [
        [1, 0], [0, 1], [1, 1], [2, 1], [3, 1], [3, 2],
        [-1, 0], [0, -1], [-1, -1], [-2, -1], [-3, -1], [-3, -2]
    ];
    for r in roots do
        if r = v then
            return true;
        fi;
    od;
    return false;
end;;

# 为 G2 场景下的 AIII 型输入构造可用根链。
# 输出的链顺序会被后续线性构造器直接消费。
BuildAIIIChainForG2 := function(ab_diagram_str, subsystem_roots, colors_opt)
    local parser, ambient, v1, v2, nc_idx, c_idx, i, nc1, nc2, vec_nc, vec_c, sum_vec, rows, rowstr, k, ch, nextch, edge_dirs, chain, ip, M;
    if not IsList(subsystem_roots) or Length(subsystem_roots) <> 2 then
        return fail;
    fi;
    parser := fail;
    ambient := "G2";
    if IsBoundGlobal("GlobalAmbientType") then
        ambient := ValueGlobal("GlobalAmbientType");
    fi;
    if ambient = "G2" then
        parser := ParseRootStringG2Simple;
    elif IsBoundGlobal("ParseRootString") then
        parser := ValueGlobal("ParseRootString");
    else
        parser := ParseRootStringG2Simple;
    fi;
    v1 := parser(subsystem_roots[1]);
    v2 := parser(subsystem_roots[2]);
    if not (IsList(v1) and IsList(v2) and Length(v1) = 2 and Length(v2) = 2) then
        return fail;
    fi;
    nc_idx := 0;
    c_idx := 0;
    if IsBoundGlobal("IsRootBlack") then
        nc1 := ValueGlobal("IsRootBlack")(v1);
        nc2 := ValueGlobal("IsRootBlack")(v2);
        if nc1 <> nc2 then
            if nc1 then
                nc_idx := 1; c_idx := 2;
            else
                nc_idx := 2; c_idx := 1;
            fi;
        fi;
    elif IsBoundGlobal("IsCompact") then
        nc1 := not ValueGlobal("IsCompact")(v1);
        nc2 := not ValueGlobal("IsCompact")(v2);
        if nc1 <> nc2 then
            if nc1 then
                nc_idx := 1; c_idx := 2;
            else
                nc_idx := 2; c_idx := 1;
            fi;
        fi;
    fi;
    if (nc_idx = 0 or c_idx = 0) and IsList(colors_opt) and Length(colors_opt) = 2 then
        for i in [1..2] do
            if colors_opt[i] = 1 then
                nc_idx := i;
            elif colors_opt[i] = 0 then
                c_idx := i;
            fi;
        od;
    fi;
    if nc_idx = 0 or c_idx = 0 then
        return fail;
    fi;
    vec_nc := parser(subsystem_roots[nc_idx]);
    vec_c := parser(subsystem_roots[c_idx]);
    sum_vec := [vec_nc[1] + vec_c[1], vec_nc[2] + vec_c[2]];
    if not IsKnownG2RootVec(sum_vec) then
        return fail;
    fi;
    if IsBoundGlobal("G2_InnerProdMatrix") then
        M := ValueGlobal("G2_InnerProdMatrix");
        ip := vec_nc[1] * M[1][1] * sum_vec[1] + vec_nc[1] * M[1][2] * sum_vec[2]
            + vec_nc[2] * M[2][1] * sum_vec[1] + vec_nc[2] * M[2][2] * sum_vec[2];
        if ip > 0 then
            sum_vec := [ -sum_vec[1], -sum_vec[2] ];
        fi;
    fi;
    rows := ParseABDiagram(ab_diagram_str);
    edge_dirs := [];
    for i in [1..Length(rows)] do
        rowstr := rows[i];
        if Length(rowstr) > 1 then
            for k in [1..Length(rowstr)-1] do
                ch := rowstr[k];
                nextch := rowstr[k+1];
                if ch = 'a' and nextch = 'b' then
                    Add(edge_dirs, 1);
                elif ch = 'b' and nextch = 'a' then
                    Add(edge_dirs, -1);
                else
                    Add(edge_dirs, 1);
                fi;
            od;
        fi;
    od;
    chain := [];
    for i in [1..Length(edge_dirs)] do
        if edge_dirs[i] = 1 then
            Add(chain, FormatRoot(sum_vec));
        else
            Add(chain, FormatRoot(vec_nc));
        fi;
    od;
    return chain;
end;;

# BDI 线性系统构造：解 A*c = b，使得对所有 e 的根 α_i，α_i(h)=2；h=∑ c_j H_{γ_j}
# 用线性代数方法构造 BDI/AIII 风格的 normal triple。
# 该函数是旧实验路径中的核心实现，仍被若干历史辅助逻辑引用。
BuildBDINormalTriple_Linear := function(ab_diagram_str, type, rank, subsystem_roots, colors_opt)
    local rows, needed, r, len, k, chain, parser, S, i, j, alpha, gamma, A, b, c, h_expr, e_expr, f_expr, allroots, v, name, seen, num, den, s, ci, sig, rot, chain2, vec1, vec2, sum_vec, root_vecs, root_names, diag_type, is_aiii, edge_dirs, rowstr, ch, nextch, FlipRoot, base_root, lhs, ok;
    rows := ParseABDiagram(ab_diagram_str);
    chain := [];
    is_aiii := (type = "su_pq" or type = "AIII" or PositionSublist(type, "su(") <> fail or PositionSublist(type, "AIII") <> fail);
    if is_aiii then
        chain := ShallowCopy(subsystem_roots);
        if Length(chain) = 0 then
            chain := ShallowCopy(subsystem_roots);
        fi;
    fi;
    if (not is_aiii) and rank = 2 and Length(subsystem_roots) = 2 and IsList(colors_opt) and Length(colors_opt) = 2 and IsBoundGlobal("ParseRootString") then
        vec1 := ValueGlobal("ParseRootString")(subsystem_roots[1]);
        vec2 := ValueGlobal("ParseRootString")(subsystem_roots[2]);
        if IsList(vec1) and IsList(vec2) and Length(vec1) = 2 and Length(vec2) = 2 then
            sum_vec := [vec1[1] + vec2[1], vec1[2] + vec2[2]];
            root_vecs := [vec1, vec2, sum_vec];
            root_names := [subsystem_roots[1], subsystem_roots[2], FormatRoot(sum_vec)];
            if IsBoundGlobal("IsRootBlack") then
                for i in [1..3] do
                    if ValueGlobal("IsRootBlack")(root_vecs[i]) then
                        Add(chain, root_names[i]);
                    fi;
                od;
            fi;
        fi;
    fi;
    if Length(chain) = 0 then chain := subsystem_roots; fi;
    # 统计需要的非紧根数
    needed := 0;
    for r in [1..Length(rows)] do
        len := Length(rows[r]);
        if len > 1 then needed := needed + (len - 1); fi;
    od;
    edge_dirs := [];
    for r in [1..Length(rows)] do
        rowstr := rows[r];
        if Length(rowstr) > 1 then
            for k in [1..Length(rowstr)-1] do
                ch := rowstr[k];
                nextch := rowstr[k+1];
                if ch = 'a' and nextch = 'b' then
                    Add(edge_dirs, 1);
                elif ch = 'b' and nextch = 'a' then
                    Add(edge_dirs, -1);
                else
                    Add(edge_dirs, 1);
                fi;
            od;
        fi;
    od;
    if is_aiii and rank = 2 and Length(subsystem_roots) = 2 then
        if IsBoundGlobal("GlobalAmbientType") and ValueGlobal("GlobalAmbientType") = "G2" then
            chain2 := BuildAIIIChainForG2(ab_diagram_str, subsystem_roots, colors_opt);
            if IsList(chain2) and Length(chain2) > 0 then
                chain := chain2;
            fi;
        fi;
    fi;
    if needed = 0 then
        return rec(h := "0", e := "0", f := "0");
    fi;
    
    sig := 0;
    for r in [1..Length(rows)] do
        sig := sig + r * Length(rows[r]);
    od;
    if (not is_aiii) and Length(chain) > 1 then
        rot := ((sig - 1) mod Length(chain)) + 1;
        chain2 := [];
        for i in [rot..Length(chain)] do Add(chain2, chain[i]); od;
        for i in [1..rot-1] do Add(chain2, chain[i]); od;
        chain := chain2;
    fi;
    # 取前 needed 个根作为 e-根链；AIII 按行内位置允许循环复用
    if Length(chain) >= needed then
        chain := chain{[1..needed]};
    elif (not is_aiii) and Length(chain) > 0 then
        needed := Length(chain);
    fi;
    # 解析器与对称矩阵
    parser := fail;
    if rank = 2 and Length(subsystem_roots) = 2 then
        parser := ParseRootStringG2Simple;
    elif IsBoundGlobal("ParseRootString") then
        parser := ValueGlobal("ParseRootString");
    else
        parser := ParseRootStringG2Simple;
    fi;
    if parser = fail then
        return rec(h := "0", e := "0", f := "0");
    fi;
    
    if rank = 2 and Length(subsystem_roots) = 2 then
        S := [[2, -3], [-3, 6]];
    else
        return rec(h := "0", e := "0", f := "0");
    fi;
    # 对根字符串整体取负，便于在线性方案中切换方向。
    FlipRoot := function(rt)
        if Length(rt) > 0 and rt[1] = '-' then
            if Length(rt) = 1 then
                return rt;
            fi;
            return rt{[2..Length(rt)]};
        fi;
        return Concatenation("-", rt);
    end;
    if (not is_aiii) and Length(chain) > 0 and needed > 0 then
        chain2 := [];
        for i in [1..needed] do
            base_root := chain[((i-1) mod Length(chain)) + 1];
            if i <= Length(edge_dirs) and edge_dirs[i] = -1 then
                Add(chain2, FlipRoot(base_root));
            else
                Add(chain2, base_root);
            fi;
        od;
        chain := chain2;
    fi;
    if IsList(chain) then
        chain := Compacted(chain);
    fi;
    if Length(chain) = 0 then
        chain := ShallowCopy(subsystem_roots);
        if IsList(chain) then
            chain := Compacted(chain);
        fi;
    fi;
    if is_aiii and Length(chain) > 0 and Length(chain) < needed then
        seen := ShallowCopy(chain);
        for i in [1..Length(subsystem_roots)] do
            if not (subsystem_roots[i] in seen) then
                Add(seen, subsystem_roots[i]);
                Add(chain, subsystem_roots[i]);
                if Length(chain) >= needed then
                    break;
                fi;
            fi;
        od;
    fi;
    if (not is_aiii) and Length(chain) > 0 and Length(chain) < needed then
        chain2 := [];
        for i in [1..needed] do
            Add(chain2, chain[((i - 1) mod Length(chain)) + 1]);
        od;
        chain := chain2;
    fi;
    if is_aiii and Length(chain) > 0 and Length(chain) < needed then
        Print("        [诊断] AIII 链长度不足，按唯一根降阶: needed=", needed, ", available=", Length(chain), "\n");
        needed := Length(chain);
    fi;
    if needed = 0 then
        return rec(h := "0", e := "0", f := "0");
    fi;
    # 构造 α_i = γ_i = 链上根（方阵尽可能接近对角）
    A := [];
    b := List([1..needed], i -> 2);
    for i in [1..needed] do
        A[i] := [];
        for j in [1..needed] do
            num := 0;
            den := 1;
            alpha := parser(chain[i]);
            gamma := parser(chain[j]);
            if IsList(alpha) and IsList(gamma) then
                if Length(alpha) = 2 and Length(gamma) = 2 and rank = 2 and Length(subsystem_roots) = 2 then
                    num := alpha[1]*(S[1][1]*gamma[1]+S[1][2]*gamma[2])
                         + alpha[2]*(S[2][1]*gamma[1]+S[2][2]*gamma[2]);
                    den := gamma[1]*(S[1][1]*gamma[1]+S[1][2]*gamma[2])
                         + gamma[2]*(S[2][1]*gamma[1]+S[2][2]*gamma[2]);
                elif Length(alpha) = 4 and Length(gamma) = 4 then
                    num := alpha[1]*(S[1][1]*gamma[1]+S[1][2]*gamma[2]+S[1][3]*gamma[3]+S[1][4]*gamma[4])
                         + alpha[2]*(S[2][1]*gamma[1]+S[2][2]*gamma[2]+S[2][3]*gamma[3]+S[2][4]*gamma[4])
                         + alpha[3]*(S[3][1]*gamma[1]+S[3][2]*gamma[2]+S[3][3]*gamma[3]+S[3][4]*gamma[4])
                         + alpha[4]*(S[4][1]*gamma[1]+S[4][2]*gamma[2]+S[4][3]*gamma[3]+S[4][4]*gamma[4]);
                    den := gamma[1]*(S[1][1]*gamma[1]+S[1][2]*gamma[2]+S[1][3]*gamma[3]+S[1][4]*gamma[4])
                         + gamma[2]*(S[2][1]*gamma[1]+S[2][2]*gamma[2]+S[2][3]*gamma[3]+S[2][4]*gamma[4])
                         + gamma[3]*(S[3][1]*gamma[1]+S[3][2]*gamma[2]+S[3][3]*gamma[3]+S[3][4]*gamma[4])
                         + gamma[4]*(S[4][1]*gamma[1]+S[4][2]*gamma[2]+S[4][3]*gamma[3]+S[4][4]*gamma[4]);
                fi;
            fi;
            if den <> 0 then
                A[i][j] := (2 * num) / den;
            else
                A[i][j] := 0;
            fi;
        od;
    od;
    # 求解（严格采用唯一线性解，不做兜底）
    c := SolveRationalLinearSystem(A, b);
    if c = fail then
        Print("        [诊断] 线性方程组无解或非唯一，返回零三元组\n");
        return rec(h := "0", e := "0", f := "0");
    fi;
    ok := true;
    for i in [1..needed] do
        lhs := 0;
        for j in [1..needed] do
            lhs := lhs + A[i][j] * c[j];
        od;
        if IsFloat(lhs) then
            if AbsFloat(lhs - 2) > 1.0e-6 then
                ok := false;
                break;
            fi;
        else
            if lhs <> 2 then
                ok := false;
                break;
            fi;
        fi;
    od;
    if not ok then
        Print("        [诊断] 线性方程组校验失败，返回零三元组\n");
        return rec(h := "0", e := "0", f := "0");
    fi;
    has_zero_coeff := false;
    for i in [1..needed] do
        if IsFloat(c[i]) then
            if AbsFloat(c[i]) <= 1.0e-6 then
                has_zero_coeff := true;
                break;
            fi;
        else
            if c[i] = 0 then
                has_zero_coeff := true;
                break;
            fi;
        fi;
    od;
    diag_type := "BDI";
    if type = "su_pq" or PositionSublist(type, "AIII") <> fail then
        diag_type := "AIII";
    fi;
    Print("        [诊断] ", diag_type, " 线性系数 c = ", c, "\n");
    # 构造表达式
    h_expr := "";
    e_expr := "";
    f_expr := "";
    for i in [1..needed] do
        if e_expr <> "" then e_expr := Concatenation(e_expr, " + "); fi;
        e_expr := Concatenation(e_expr, "E_{[", chain[i], "]}");
    od;
    for i in [1..needed] do
        if (IsFloat(c[i]) and AbsFloat(c[i]) <= 1.0e-6) or ((not IsFloat(c[i])) and c[i] = 0) then
            continue;
        fi;
        if h_expr <> "" then h_expr := Concatenation(h_expr, " + "); fi;
        if c[i] = 1 then
            h_expr := Concatenation(h_expr, "H_{", chain[i], "}");
        else
            h_expr := Concatenation(h_expr, String(c[i]), "H_{", chain[i], "}");
        fi;
        if f_expr <> "" then f_expr := Concatenation(f_expr, " + "); fi;
        if c[i] = 1 then
            f_expr := Concatenation(f_expr, "F_{[", chain[i], "]}");
        else
            f_expr := Concatenation(f_expr, String(c[i]), "F_{[", chain[i], "]}");
        fi;
    od;
    if h_expr = "" then h_expr := "0"; fi;
    if e_expr = "" then e_expr := "0"; fi;
    if f_expr = "" then f_expr := "0"; fi;
    Print("        [LINEAR-CANDIDATE] h = ", h_expr, "\n");
    Print("        [LINEAR-CANDIDATE] e = ", e_expr, "\n");
    Print("        [LINEAR-CANDIDATE] f = ", f_expr, "\n");
    return rec(h := h_expr, e := e_expr, f := f_expr);
end;;

# 模块 4 对外暴露的统一打印接口。
# 当前实现直接转发到统一主入口 PrintSL2TripleUnified。
PrintSL2Triple := function(ab_diagram_str, type, rank, subsystem_roots, colors_opt)
    PrintSL2TripleUnified(ab_diagram_str, type, rank, subsystem_roots, colors_opt);
end;;

# 解析形如 "sqrt(6)E_{[a1]}" 或 "2H_{a1}" 的项，返回记录 {kind, root, coeff}
# 解析表达式中的单个 H/E/F 项。
# 返回项种类、根名与系数，供校验与重整逻辑复用。
ParseTerm := function(term)
    local i, kind, root, coeff_str, coeff, idx_H, idx_E, idx_F, s, start_idx, end_idx, ParseRational, parsed_val;
    # 解析整数或分数字符串为有理数。
    ParseRational := function(str)
        local slash, num_str, den_str, num, den, val;
        slash := PositionSublist(str, "/");
        if slash <> fail then
            num_str := str{[1..slash-1]};
            den_str := str{[slash+1..Length(str)]};
            num := Int(String(num_str));
            den := Int(String(den_str));
            if num = fail or den = fail or den = 0 then
                return fail;
            fi;
            return num/den;
        fi;
        val := Int(String(str));
        if val <> fail then
            return val;
        fi;
        if PositionSublist(str, ".") <> fail or PositionSublist(str, "e") <> fail or PositionSublist(str, "E") <> fail then
            val := EvalString(str);
            if IsRat(val) or IsFloat(val) then
                return val;
            fi;
        fi;
        return fail;
    end;
    s := term;
    kind := "";
    root := "";
    coeff := 1;
    idx_H := PositionSublist(s, "H_{");
    idx_E := PositionSublist(s, "E_{");
    idx_F := PositionSublist(s, "F_{");
    if idx_H <> fail then kind := "H";
    elif idx_E <> fail then kind := "E";
    elif idx_F <> fail then kind := "F";
    fi;
    if kind = "" then
        return rec(kind := "", root := "", coeff := 0);
    fi;
    if PositionSublist(s, "1/sqrt(") = 1 then
        start_idx := 8;
        end_idx := PositionSublist(s, ")");
        if end_idx = fail then
            coeff := 1;
        else
            end_idx := end_idx - 1;
            coeff_str := s{[start_idx..end_idx]};
            coeff := ParseRational(coeff_str);
            if coeff = fail then
                coeff := 1;
            elif IsInt(coeff) then
                coeff := 1 / Sqrt(coeff);
            else
                coeff := 1 / Sqrt(Float(coeff));
            fi;
        fi;
    elif PositionSublist(s, "sqrt(") = 1 then
        start_idx := 6;
        end_idx := PositionSublist(s, ")");
        if end_idx = fail then
            coeff := 1;
        else
            end_idx := end_idx - 1;
            coeff_str := s{[start_idx..end_idx]};
            coeff := ParseRational(coeff_str);
            if coeff = fail then
                coeff := 1;
            elif IsInt(coeff) then
                coeff := Sqrt(coeff);
            else
                coeff := Sqrt(Float(coeff));
            fi;
        fi;
    else
        # 可能以数字系数开头，例如 "2H_{...}"
        i := 1;
        while i <= Length(s) and s[i] in "0123456789-/.eE+" do
            i := i + 1;
        od;
        if i > 1 then
            coeff_str := s{[1..i-1]};
            coeff := ParseRational(coeff_str);
            if coeff = fail then coeff := 1; fi;
        else
            coeff := 1;
        fi;
    fi;
    start_idx := PositionSublist(s, "{");
    end_idx := PositionSublist(s, "}");
    if start_idx = fail or end_idx = fail or end_idx <= start_idx then
        return rec(kind := "", root := "", coeff := 0);
    fi;
    root := s{[start_idx+1..end_idx-1]};
    if Length(root) >= 2 and root[1] = '[' and root[Length(root)] = ']' then
        root := root{[2..Length(root)-1]};
    fi;
    return rec(kind := kind, root := root, coeff := coeff);
end;;

# G2 根向量格式化（简单根基 a1, a2）
# 将 G2 根向量格式化为字符串。
# 输出形如 "a1+a2"、"-3a1-a2"，供日志与表达式拼接使用。
FormatRoot := function(v)
    local s, i, coeff, term;
    s := "";
    for i in [1..2] do
        coeff := v[i];
        if coeff <> 0 then
            if s <> "" and coeff > 0 then s := Concatenation(s, "+"); fi;
            if coeff = 1 then
                term := Concatenation("a", String(i));
            elif coeff = -1 then
                term := Concatenation("-a", String(i));
            else
                term := Concatenation(String(coeff), "a", String(i));
            fi;
            s := Concatenation(s, term);
        fi;
    od;
    if s = "" then s := "0"; fi;
    return s;
end;;

# 将形如 "sqrt(6)E_{[a1]} + E_{[a1+a2]}" 的字符串解析为 root->系数之和的字典（以列表对实现）
# 再次定义表达式系数提取器，供后段校验逻辑独立使用。
# 它把表达式拆成“根名 → 系数”的映射，便于比较左右两边。
CollectCoeffs := function(expr, expected_kind)
    local parts, p, term, data, acc, found, pair, i, Trim, rest, pos;
    # 简单去除首尾空格
    # 去掉项字符串的首尾空白。
    Trim := function(s)
        local l, r;
        l := 1; r := Length(s);
        while l <= r and s[l] = ' ' do l := l + 1; od;
        while r >= l and s[r] = ' ' do r := r - 1; od;
        if r < l then return ""; fi;
        return s{[l..r]};
    end;
    acc := [];
    if expr = "0" then return acc; fi;
    parts := [];
    rest := expr;
    while true do
        pos := PositionSublist(rest, " + ");
        if pos = fail then
            Add(parts, rest);
            break;
        fi;
        Add(parts, rest{[1..pos-1]});
        if pos + 3 <= Length(rest) then
            rest := rest{[pos+3..Length(rest)]};
        else
            rest := "";
        fi;
    od;
    for p in parts do
        term := Trim(p);
        data := ParseTerm(term);
        if data.kind = expected_kind then
            found := false;
            for pair in acc do
                if pair[1] = data.root then
                    pair[2] := pair[2] + data.coeff;
                    found := true;
                    break;
                fi;
            od;
            if not found then
                Add(acc, [ data.root, data.coeff ]);
            fi;
        fi;
    od;
    return acc;
end;;

# 校验函数：打印 [e,f]=h 是否成立；以及简单的 [h,e]=2e、[h,f]=-2f 系数检查
# 验证 sl2-triple 的 EFH 正常性条件。
# 包括 [e,f]=h、[h,e]=2e、[h,f]=-2f 以及表达式层面的结构检查。
VerifyEFH_NormalConditions := function(h_str, e_str, f_str)
    local h_map, e_terms, f_terms, i, n, accum, pair, root, ch, sumv, ok, pr, he_ok, hf_ok, parser, alpha_h, cond, logf, h_terms, eval_alpha, ParseTermEFH, NormalizeExprForParse, SplitExprTerms, h_norm, e_norm, f_norm, h_coeffs, MirrorEFHText;
    if IsBound(GlobalOutFile) then
        logf := GlobalOutFile;
    elif IsBound(GlobalOutputDir) then
        logf := Concatenation(GlobalOutputDir, "/G2_pipeline_output.txt");
    else
        logf := "G2_pipeline_output.txt";
    fi;
    MirrorEFHText := function(text)
        if IsBoundGlobal("AppendRetainedMirrorText") then
            ValueGlobal("AppendRetainedMirrorText")(text);
        fi;
    end;
    # 预清洗表达式，便于后续逐项解析。
    NormalizeExprForParse := function(expr)
        local s;
        s := Filtered(expr, c -> c <> '\n' and c <> '\r' and c <> '\t' and c <> '\\');
        while PositionSublist(s, "  ") <> fail do
            s := ReplacedString(s, "  ", " ");
        od;
        return s;
    end;
    h_norm := NormalizeExprForParse(h_str);
    e_norm := NormalizeExprForParse(e_str);
    f_norm := NormalizeExprForParse(f_str);
    # 把表达式按顶层加号拆成项列表。
    SplitExprTerms := function(expr)
        local terms, rest, pos;
        terms := [];
        rest := expr;
        while true do
            pos := PositionSublist(rest, " + ");
            if pos = fail then
                Add(terms, rest);
                break;
            fi;
            Add(terms, rest{[1..pos-1]});
            if pos + 3 <= Length(rest) then
                rest := rest{[pos+3..Length(rest)]};
            else
                rest := "";
            fi;
        od;
        return terms;
    end;
    # 有些构造会在同一根上出现多次项，此时不能先求和再相乘（会引入交叉项）
    # 采用按序配对的方法：对第 i 项计算 e_i * f_i，并按根名累加
    # 解析 EFH 校验阶段使用的单项字符串。
    ParseTermEFH := function(term)
        local s, kind, root, coeff, idx_H, idx_E, idx_F, start_idx, end_idx, i2, coeff_str, eff_coeff, ParseRational, slash, num_str1, num_str2, num1, num2, num_val, parsed_val;
        # 解析 EFH 局部项中的有理系数。
        ParseRational := function(str)
            slash := PositionSublist(str, "/");
            if slash <> fail then
                num_str1 := str{[1..slash-1]};
                num_str2 := str{[slash+1..Length(str)]};
                num1 := Int(String(num_str1));
                num2 := Int(String(num_str2));
                if num1 = fail or num2 = fail or num2 = 0 then
                    return fail;
                fi;
                return num1/num2;
            fi;
            parsed_val := Int(String(str));
            if parsed_val <> fail then
                return parsed_val;
            fi;
            if PositionSublist(str, ".") <> fail or PositionSublist(str, "e") <> fail or PositionSublist(str, "E") <> fail then
                parsed_val := EvalString(str);
                if IsRat(parsed_val) or IsFloat(parsed_val) then
                    return parsed_val;
                fi;
            fi;
            return fail;
        end;
        s := term;
        kind := ""; root := ""; coeff := 1; eff_coeff := fail;
        idx_H := PositionSublist(s, "H_{");
        idx_E := PositionSublist(s, "E_{");
        idx_F := PositionSublist(s, "F_{");
        if idx_H <> fail then kind := "H";
        elif idx_E <> fail then kind := "E";
        elif idx_F <> fail then kind := "F";
        fi;
        if kind = "" then
            return rec(kind := "", root := "", coeff := 0, eff_coeff := fail);
        fi;
        if PositionSublist(s, "1/sqrt(") = 1 then
            start_idx := 8;
            end_idx := PositionSublist(s, ")");
            if end_idx <> fail then
                end_idx := end_idx - 1;
                coeff_str := s{[start_idx..end_idx]};
                num_val := ParseRational(coeff_str);
                if num_val <> fail and num_val <> 0 then
                    eff_coeff := 1 / num_val;
                    coeff := 1 / Sqrt(Float(num_val));
                fi;
            fi;
        elif PositionSublist(s, "sqrt(") = 1 then
            start_idx := 6;
            end_idx := PositionSublist(s, ")");
            if end_idx <> fail then
                end_idx := end_idx - 1;
                coeff_str := s{[start_idx..end_idx]};
                num_val := ParseRational(coeff_str);
                if num_val <> fail then
                    eff_coeff := num_val;
                    coeff := Sqrt(Float(num_val));
                fi;
            fi;
        else
            i2 := 1;
            while i2 <= Length(s) and s[i2] in "0123456789-/.eE+" do
                i2 := i2 + 1;
            od;
            if i2 > 1 then
                coeff_str := s{[1..i2-1]};
                num_val := ParseRational(coeff_str);
                if num_val = fail then num_val := 1; fi;
                coeff := num_val;
                eff_coeff := num_val * num_val;
            else
                coeff := 1;
                eff_coeff := 1;
            fi;
        fi;
        start_idx := PositionSublist(s, "{");
        end_idx := PositionSublist(s, "}");
        if start_idx = fail or end_idx = fail or end_idx <= start_idx then
            return rec(kind := "", root := "", coeff := 0, eff_coeff := fail);
        fi;
        root := s{[start_idx+1..end_idx-1]};
        if Length(root) >= 2 and root[1] = '[' and root[Length(root)] = ']' then
            root := root{[2..Length(root)-1]};
        fi;
        return rec(kind := kind, root := root, coeff := coeff, eff_coeff := eff_coeff);
    end;
    e_terms := []; f_terms := [];
    if e_norm <> "0" then
        for pr in SplitExprTerms(e_norm) do
            pair := ParseTermEFH(pr);
            if pair.kind = "E" then Add(e_terms, [pair.root, pair.coeff, pair.eff_coeff]); fi;
        od;
    fi;
    if f_norm <> "0" then
        for pr in SplitExprTerms(f_norm) do
            pair := ParseTermEFH(pr);
            if pair.kind = "F" then Add(f_terms, [pair.root, pair.coeff, pair.eff_coeff]); fi;
        od;
    fi;
    # 目标 h 的聚合（按根名求和）
    h_map := CollectCoeffs(h_norm, "H");
    # 逐项相乘累加
    accum := [];
    n := Minimum(Length(e_terms), Length(f_terms));
    for i in [1..n] do
        root := e_terms[i][1];
        if root <> f_terms[i][1] then
            # 根名不一致，无法匹配，直接失败
            Print("        [诊断] 第 ", i, " 项根不匹配: E=", e_terms[i][1], " F=", f_terms[i][1], "\n");
            MirrorEFHText(Concatenation("        [诊断] 第 ", String(i), " 项根不匹配: E=", e_terms[i][1], " F=", f_terms[i][1], "\n"));
            AppendTo(logf, Concatenation("        [诊断] 第 ", String(i), " 项根不匹配: E=", e_terms[i][1], " F=", f_terms[i][1], "\n"));
            Print("EFH_CHECK: FAIL\n");
            MirrorEFHText("EFH_CHECK: FAIL\n");
            AppendTo(logf, "EFH_CHECK: FAIL\n");
            return;
        fi;
        if e_terms[i][3] <> fail and f_terms[i][3] <> fail and e_terms[i][3] = f_terms[i][3] then
            sumv := e_terms[i][3];
        else
            sumv := e_terms[i][2] * f_terms[i][2];
        fi;
        # 累加到对应根
        ok := false;
        for pair in accum do
            if pair[1] = root then
                pair[2] := pair[2] + sumv;
                ok := true;
                break;
            fi;
        od;
        if not ok then
            Add(accum, [root, sumv]);
        fi;
    od;
    # 对比 accum 与 h_map
    ok := true;
    for pair in accum do
        root := pair[1];
        sumv := pair[2];
        ch := 0;
        for pr in h_map do
            if pr[1] = root then
                ch := pr[2];
                break;
            fi;
        od;
        # 采用容差比较，避免浮点误差
        if IsFloat(sumv) or IsFloat(ch) then
            if AbsFloat(Float(sumv) - Float(ch)) > 1.0e-6 then ok := false; fi;
        else
            if sumv <> ch then ok := false; fi;
        fi;
    od;
    if ok then
        Print("        校验: [e,f]=h: PASS\n");
        MirrorEFHText("        校验: [e,f]=h: PASS\n");
        AppendTo(logf, "        校验: [e,f]=h: PASS\n");
        Print("EFH_CHECK: PASS\n");
        MirrorEFHText("EFH_CHECK: PASS\n");
        AppendTo(logf, "EFH_CHECK: PASS\n");
    else
        Print("        校验: [e,f]=h: FAIL\n");
        MirrorEFHText("        校验: [e,f]=h: FAIL\n");
        AppendTo(logf, "        校验: [e,f]=h: FAIL\n");
        Print("EFH_CHECK: FAIL\n");
        MirrorEFHText("EFH_CHECK: FAIL\n");
        AppendTo(logf, "EFH_CHECK: FAIL\n");
    fi;
    he_ok := true;
    hf_ok := true;
    if IsBoundGlobal("ParseRootString") or IsBoundGlobal("ParseRootStringG2Simple") then
        parser := ParseRootStringG2Simple;
        if IsBoundGlobal("ParseRootString") then
            parser := ValueGlobal("ParseRootString");
        fi;
        h_terms := CollectCoeffs(h_norm, "H");
        h_coeffs := fail;
        if IsBoundGlobal("ParseHString") then
            h_coeffs := ValueGlobal("ParseHString")(h_norm);
        fi;
        # 计算单根 alpha 在当前 h 下的权值 alpha(h)。
        eval_alpha := function(alpha_root_str)
            local v, total, pair, g, num, den;
            v := parser(alpha_root_str);
            if h_coeffs <> fail and IsList(v) and Length(v) >= 2 then
                if IsBoundGlobal("CalculateWeight") then
                    return ValueGlobal("CalculateWeight")(v, h_coeffs);
                fi;
                return v[1] * (2*h_coeffs[1] - h_coeffs[2]) + v[2] * (-3*h_coeffs[1] + 2*h_coeffs[2]);
            fi;
            total := 0;
            for pair in h_terms do
                g := parser(pair[1]);
                if IsList(v) and IsList(g) and Length(v) >= 2 and Length(g) >= 2 then
                    num := v[1]*(2*g[1]-3*g[2]) + v[2]*(-3*g[1]+6*g[2]);
                    den := g[1]*(2*g[1]-3*g[2]) + g[2]*(-3*g[1]+6*g[2]);
                else
                    num := 0;
                    den := 0;
                fi;
                if den <> 0 then
                    total := total + pair[2] * (2 * num / den);
                fi;
            od;
            return total;
        end;
        for i in [1..Length(e_terms)] do
            alpha_h := eval_alpha(e_terms[i][1]);
            if IsFloat(alpha_h) then
                cond := AbsFloat(alpha_h - 2) <= 1.0e-6;
            else
                cond := (alpha_h = 2);
            fi;
            if cond then
                Print("        [诊断] E 根 ", e_terms[i][1], " 的 α(h)=", alpha_h, " : PASS\n");
                MirrorEFHText(Concatenation("        [诊断] E 根 ", e_terms[i][1], " 的 α(h)=", String(alpha_h), " : PASS\n"));
                AppendTo(logf, Concatenation("        [诊断] E 根 ", e_terms[i][1], " 的 α(h)=", String(alpha_h), " : PASS\n"));
            else
                Print("        [诊断] E 根 ", e_terms[i][1], " 的 α(h)=", alpha_h, " : FAIL\n");
                MirrorEFHText(Concatenation("        [诊断] E 根 ", e_terms[i][1], " 的 α(h)=", String(alpha_h), " : FAIL\n"));
                AppendTo(logf, Concatenation("        [诊断] E 根 ", e_terms[i][1], " 的 α(h)=", String(alpha_h), " : FAIL\n"));
                he_ok := false;
                break;
            fi;
        od;
        for i in [1..Length(f_terms)] do
            alpha_h := -eval_alpha(f_terms[i][1]);
            if IsFloat(alpha_h) then
                cond := AbsFloat(alpha_h + 2) <= 1.0e-6;
            else
                cond := (alpha_h = -2);
            fi;
            if cond then
                Print("        [诊断] F 根 ", f_terms[i][1], " 的 α(h)=", -alpha_h, " => 检查 -α(h)=-2: PASS\n");
                MirrorEFHText(Concatenation("        [诊断] F 根 ", f_terms[i][1], " 的 α(h)=", String(-alpha_h), " => 检查 -α(h)=-2: PASS\n"));
                AppendTo(logf, Concatenation("        [诊断] F 根 ", f_terms[i][1], " 的 α(h)=", String(-alpha_h), " => 检查 -α(h)=-2: PASS\n"));
            else
                Print("        [诊断] F 根 ", f_terms[i][1], " 的 α(h)=", -alpha_h, " => 检查 -α(h)=-2: FAIL\n");
                MirrorEFHText(Concatenation("        [诊断] F 根 ", f_terms[i][1], " 的 α(h)=", String(-alpha_h), " => 检查 -α(h)=-2: FAIL\n"));
                AppendTo(logf, Concatenation("        [诊断] F 根 ", f_terms[i][1], " 的 α(h)=", String(-alpha_h), " => 检查 -α(h)=-2: FAIL\n"));
                hf_ok := false;
                break;
            fi;
        od;
        if he_ok then
            Print("        校验: [h,e]=2e: PASS\n");
            MirrorEFHText("        校验: [h,e]=2e: PASS\n");
            AppendTo(logf, "        校验: [h,e]=2e: PASS\n");
        else
            Print("        校验: [h,e]=2e: FAIL\n");
            MirrorEFHText("        校验: [h,e]=2e: FAIL\n");
            AppendTo(logf, "        校验: [h,e]=2e: FAIL\n");
        fi;
        if hf_ok then
            Print("        校验: [h,f]=-2f: PASS\n");
            MirrorEFHText("        校验: [h,f]=-2f: PASS\n");
            AppendTo(logf, "        校验: [h,f]=-2f: PASS\n");
        else
            Print("        校验: [h,f]=-2f: FAIL\n");
            MirrorEFHText("        校验: [h,f]=-2f: FAIL\n");
            AppendTo(logf, "        校验: [h,f]=-2f: FAIL\n");
        fi;
        if ok and he_ok and hf_ok then
            return true;
        else
            return false;
        fi;
    else
        if ok then
            return true;
        else
            return false;
        fi;
    fi;
end;;

# 验证 G2 子系统上的 P/K 归属条件。
# 要求 E、F 项位于非紧部分，H 项与子系统结构保持一致。
VerifyPK_Subsystem_G2 := function(h_str, e_str, f_str, subsystem_roots, colors_opt)
    local i, root_name, nc_flag, ok, ok_e, ok_f, ok_h, ParseRootSimple, IsNoncompactRoot, CheckTerms, p, data, ComputeSubsystemCoeffs, black_sum;
    ParseRootSimple := ParseRootStringG2Simple;
    # 计算根在当前子系统基下的系数表示。
    ComputeSubsystemCoeffs := function(root_str)
        local v, r1, r2, det, c1, c2;
        if not IsList(subsystem_roots) or Length(subsystem_roots) < 2 then
            return fail;
        fi;
        if IsBoundGlobal("ParseRootString") then
            v := ValueGlobal("ParseRootString")(root_str);
            r1 := ValueGlobal("ParseRootString")(subsystem_roots[1]);
            r2 := ValueGlobal("ParseRootString")(subsystem_roots[2]);
        else
            v := ParseRootSimple(root_str);
            r1 := ParseRootSimple(subsystem_roots[1]);
            r2 := ParseRootSimple(subsystem_roots[2]);
        fi;
        det := r1[1] * r2[2] - r1[2] * r2[1];
        if det = 0 then
            return fail;
        fi;
        c1 := (v[1] * r2[2] - v[2] * r2[1]) / det;
        c2 := (r1[1] * v[2] - r1[2] * v[1]) / det;
        if IsInt(c1) and IsInt(c2) then
            return [c1, c2];
        fi;
        return fail;
    end;
    # 判断字符串形式的根是否属于非紧部分。
    IsNoncompactRoot := function(root_str)
        local vec, is_black, is_compact;
        if IsBoundGlobal("ParseRootString") then
            vec := ValueGlobal("ParseRootString")(root_str);
        else
            vec := ParseRootSimple(root_str);
        fi;
        if IsBoundGlobal("IsRootBlack") then
            is_black := ValueGlobal("IsRootBlack")(vec);
            return is_black;
        elif IsBoundGlobal("IsCompact") then
            is_compact := ValueGlobal("IsCompact")(vec);
            return not is_compact;
        fi;
        if IsList(vec) and Length(vec) >= 2 then
            if vec[2] < 0 then vec[2] := -vec[2]; fi;
            return (vec[2] mod 2) <> 0;
        fi;
        return fail;
    end;
    # 检查表达式各项是否都满足 P/K 归属要求。
    CheckTerms := function(expr, expect_noncompact)
        local parts, term, nc, NormalizeRootString, root_clean, SplitByPlusOutside, buf, ch, Trim, depth;
        # 标准化根字符串，便于比较。
        NormalizeRootString := function(s)
            return Filtered(s, c -> c <> ' ' and c <> '[' and c <> ']');
        end;
        # 去除局部字符串首尾空白。
        Trim := function(s)
            local l, r;
            l := 1; r := Length(s);
            while l <= r and s[l] = ' ' do l := l + 1; od;
            while r >= l and s[r] = ' ' do r := r - 1; od;
            if r < l then return ""; fi;
            return s{[l..r]};
        end;
        # 在不破坏括号结构的前提下按加号拆项。
        SplitByPlusOutside := function(s)
            local out, tmp;
            out := [];
            tmp := "";
            depth := 0;
            for ch in s do
                if ch = '[' or ch = '{' then
                    depth := depth + 1;
                    tmp := Concatenation(tmp, [ch]);
                elif ch = ']' or ch = '}' then
                    depth := depth - 1;
                    tmp := Concatenation(tmp, [ch]);
                elif ch = '+' and depth = 0 then
                    Add(out, Trim(tmp));
                    tmp := "";
                else
                    tmp := Concatenation(tmp, [ch]);
                fi;
            od;
            Add(out, Trim(tmp));
            return out;
        end;
        if expr = "0" then return true; fi;
        parts := SplitByPlusOutside(expr);
        for term in parts do
            if term = "" then
                Print("        [诊断] P/K 解析失败: term=\n");
                return false;
            fi;
            data := ParseTerm(term);
            if data.kind = "" then
                Print("        [诊断] P/K 解析失败: term=", term, "\n");
                return false;
            fi;
            root_clean := NormalizeRootString(data.root);
            nc := IsNoncompactRoot(root_clean);
            if nc = fail then
                Print("        [诊断] P/K 根无法判定: root=", root_clean, "\n");
                return false;
            fi;
            if expect_noncompact and nc = false then
                Print("        [诊断] P/K 根应非紧但为紧: root=", root_clean, "\n");
                return false;
            fi;
            if (not expect_noncompact) and nc = true then
                Print("        [诊断] P/K 根应紧但为非紧: root=", root_clean, "\n");
                return false;
            fi;
        od;
        return true;
    end;
    ok_e := CheckTerms(e_str, true);
    ok_f := CheckTerms(f_str, true);
    if ok_e and ok_f then
        return true;
    else
        return false;
    fi;
end;;
