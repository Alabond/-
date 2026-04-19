#############################################################################
# 文件: SL2_Triple_Builder.g
# 描述: 基于 theta-orbit 翻译与校验 sl2-triple (h, e, f)
# 说明: 该文件较大，按“解析/规范化/构造/过滤/验证/显示”六类组织函数。
# 核心入口函数用途:
#   - PrintSL2TripleUnified:
#       模块4主入口。根据类型与实形式路线，统一输出可验证的 h,e,f。
#   - BuildCompositeNormalTriple_FromComponents:
#       对可分解子系统按分量分别构造 triple，再合成总体 triple。
#   - FilterTriplesByWKConjugacyFixingDeltaJF4 /
#     FilterCompositeTripleFamiliesByKOrbitLabelF4:
#       对候选 triple 做 W_k(Δ_J 固定) 与 KOrbitLabel 层面的去重过滤。
#   - VerifyThetaSelectedTripleOnLieLevelF4 /
#     RefitCompositeTripleF4_PostCheck:
#       将候选 triple 拉回 Lie 层做一致性复核，必要时回代修正。
#   - ComputeThetaOrbitTriplesForModule4 / BuildDirectThetaTripleForF4:
#       计算 θ-轨道候选并尝试直接生成根级 triple。
#############################################################################

# =============================================================================
# 基础辅助函数（先定义，避免未绑定）
# =============================================================================

# 浮点绝对值
AbsFloat := function(x)
    local y;
    if IsFloat(x) then
        y := x;
    else
        y := Float(x);
    fi;
    if y >= 0.0 then return y; else return -y; fi;
end;;

GlobalNormalTripleTotalChecks := 0;
GlobalNormalTriplePassChecks := 0;
GlobalF4RootSet_SimpleCoeffs := fail;
GlobalCurrentThetaOrbitData := fail;

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

# 统一入口：根据需要规范化根，再复用通用构造
F4InnerProductBySym := function(v1, v2)
    local S, i, j, acc;
    if not IsBoundGlobal("F4_SymMatrix") then
        return fail;
    fi;
    S := ValueGlobal("F4_SymMatrix")();
    acc := 0;
    for i in [1..Length(v1)] do
        for j in [1..Length(v2)] do
            acc := acc + v1[i] * S[i][j] * v2[j];
        od;
    od;
    return acc;
end;;

ComputeCartanMatrixForF4RootList := function(root_vecs)
    local A, i, j, den;
    A := [];
    for i in [1..Length(root_vecs)] do
        A[i] := [];
        den := F4InnerProductBySym(root_vecs[i], root_vecs[i]);
        if den = 0 then
            return fail;
        fi;
        for j in [1..Length(root_vecs)] do
            A[i][j] := 2 * F4InnerProductBySym(root_vecs[i], root_vecs[j]) / den;
        od;
    od;
    return A;
end;;

IsIntegralSpanMemberF4 := function(basis_vecs, target_vec)
    local gram, rhs, i, j, coeffs, rebuilt, solve_linear;
    if not (IsList(basis_vecs) and Length(basis_vecs) > 0 and IsList(target_vec)) then
        return false;
    fi;
    gram := [];
    rhs := [];
    for i in [1..Length(basis_vecs)] do
        gram[i] := [];
        rhs[i] := F4InnerProductBySym(target_vec, basis_vecs[i]);
        for j in [1..Length(basis_vecs)] do
            gram[i][j] := F4InnerProductBySym(basis_vecs[j], basis_vecs[i]);
        od;
    od;
    if not IsBoundGlobal("SolveRationalLinearSystem") then
        return false;
    fi;
    solve_linear := ValueGlobal("SolveRationalLinearSystem");
    coeffs := solve_linear(gram, rhs);
    if coeffs = fail then
        return false;
    fi;
    for i in [1..Length(coeffs)] do
        if coeffs[i] <> Int(coeffs[i]) then
            return false;
        fi;
    od;
    rebuilt := List([1..Length(target_vec)], k -> 0);
    for i in [1..Length(basis_vecs)] do
        rebuilt := List([1..Length(rebuilt)], k -> rebuilt[k] + Int(coeffs[i]) * basis_vecs[i][k]);
    od;
    return rebuilt = target_vec;
end;;

VerifyC3PositiveClosureInF4 := function(simple_vecs)
    local pos_coeffs, coeffs, vec;
    pos_coeffs := [
        [1,0,0], [0,1,0], [0,0,1], [1,1,0], [0,1,1],
        [1,1,1], [0,2,1], [1,2,1], [2,2,1]
    ];
    for coeffs in pos_coeffs do
        vec := List([1..Length(simple_vecs[1])], k -> 0);
        vec := List([1..Length(vec)], k ->
            coeffs[1] * simple_vecs[1][k] +
            coeffs[2] * simple_vecs[2][k] +
            coeffs[3] * simple_vecs[3][k]
        );
        if not ValueGlobal("IsF4RootVec")(vec) then
            return false;
        fi;
    od;
    return true;
end;;

FindC3SimpleRootsInSpanF4 := function(root_vecs)
    local all_roots, span_roots, v, i, j, k, cand_vecs, A, out, target_cartans;
    if not (IsList(root_vecs) and Length(root_vecs) = 3 and IsBoundGlobal("F4_GenerateRootSet") and IsBoundGlobal("FormatRootF4")) then
        return fail;
    fi;
    all_roots := ValueGlobal("F4_GenerateRootSet")();
    target_cartans := [
        [[2,-1,0],[-1,2,-1],[0,-2,2]],
        [[2,-1,0],[-1,2,-2],[0,-1,2]]
    ];
    span_roots := [];
    for v in all_roots do
        if IsIntegralSpanMemberF4(root_vecs, v) then
            Add(span_roots, v);
        fi;
    od;
    for i in [1..Length(span_roots)] do
        for j in [1..Length(span_roots)] do
            for k in [1..Length(span_roots)] do
                if i <> j and i <> k and j <> k then
                    cand_vecs := [span_roots[i], span_roots[j], span_roots[k]];
                    A := ComputeCartanMatrixForF4RootList(cand_vecs);
                    if A in target_cartans and VerifyC3PositiveClosureInF4(cand_vecs) then
                        out := [];
                        Add(out, ValueGlobal("FormatRootF4")(cand_vecs[1]));
                        Add(out, ValueGlobal("FormatRootF4")(cand_vecs[2]));
                        Add(out, ValueGlobal("FormatRootF4")(cand_vecs[3]));
                        return out;
                    fi;
                fi;
            od;
        od;
    od;
    return fail;
end;;

NormalizeC3RootsForThetaF4 := function(roots)
    local parser, perms, signs, root_vecs, perm, sign_choice, cand_vecs, i, A, target, out, span_choice;
    if not (IsList(roots) and Length(roots) = 3) then
        return fail;
    fi;
    if IsBoundGlobal("ParseRootStringFlexible") then
        parser := ValueGlobal("ParseRootStringFlexible");
    elif IsBoundGlobal("ParseRootString") then
        parser := ValueGlobal("ParseRootString");
    else
        return fail;
    fi;
    if not IsBoundGlobal("FormatRootF4") or not IsBoundGlobal("IsF4RootVec") then
        return fail;
    fi;
    root_vecs := List(roots, r -> parser(r));
    if fail in root_vecs then
        return fail;
    fi;
    perms := [
        [1,2,3], [1,3,2], [2,1,3],
        [2,3,1], [3,1,2], [3,2,1]
    ];
    signs := [
        [1,1,1], [1,1,-1], [1,-1,1], [1,-1,-1],
        [-1,1,1], [-1,1,-1], [-1,-1,1], [-1,-1,-1]
    ];
    target := [[2,-1,0],[-1,2,-1],[0,-2,2]];
    for perm in perms do
        for sign_choice in signs do
            cand_vecs := [];
            for i in [1..3] do
                Add(cand_vecs, List(root_vecs[perm[i]], x -> sign_choice[i] * x));
            od;
            if ForAll(cand_vecs, v -> ValueGlobal("IsF4RootVec")(v)) then
                A := ComputeCartanMatrixForF4RootList(cand_vecs);
                if A = target then
                    out := [];
                    for i in [1..3] do
                        Add(out, ValueGlobal("FormatRootF4")(cand_vecs[i]));
                    od;
                    return out;
                fi;
            fi;
        od;
    od;
    span_choice := FindC3SimpleRootsInSpanF4(root_vecs);
    if span_choice <> fail then
        return span_choice;
    fi;
    return fail;
end;;

NormalizeRootsForNormalTriple := function(type, rank, roots, colors_opt)
    local i, has_compact, has_noncompact, compact_idx, non_compact_idx,
          vec_C, vec_NC, valid_g2_roots, is_root, r1_vec, r2_vec, r1_str, r2_str, sum_vec,
          M, len1, len2, ip, ambient, ParseRootSimple, ParseRoot, v1, v2, nc1, nc2, normalized_c3;
    
    if not IsList(roots) or Length(roots) = 0 then
        return roots;
    fi;

    if rank = 3 and Length(roots) = 3 and (PositionSublist(type, "sp") <> fail or PositionSublist(type, "C") = 1) then
        normalized_c3 := NormalizeC3RootsForThetaF4(roots);
        if normalized_c3 <> fail then
            return normalized_c3;
        fi;
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
            ambient := "F4";
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

CollectCoeffs := function(expr, expected_kind)
    local parts, p, term, data, acc, found, pair, Trim;
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
    parts := SplitString(expr, " + ");
    for p in parts do
        term := Trim(p);
        data := ValueGlobal("ParseTerm")(term);
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

RefitCompositeTripleForF4 := function(triple)
    local parser, S, e_terms, roots, term, pair, seen, v, n, A, b, i, j, vi, vj, num, den, c, h_expr, e_expr, f_expr, coeff, SplitByPlus, ok, lhs, e_roots, SolveUniqueSystem, BuildSystemAndSolve, ValidateAgainstERoots, c_try, ParseTermLocal, best_roots, best_try, chosen_roots, DFSSubsetSearch, subset_size, validated_ok;
    SplitByPlus := function(expr)
        local out, cur, p;
        out := [];
        cur := expr;
        while true do
            p := PositionSublist(cur, " + ");
            if p = fail then
                Add(out, cur);
                break;
            fi;
            if p > 1 then
                Add(out, cur{[1..p-1]});
            else
                Add(out, "");
            fi;
            if p + 3 <= Length(cur) then
                cur := cur{[p+3..Length(cur)]};
            else
                cur := "";
                Add(out, cur);
                break;
            fi;
        od;
        return out;
    end;
    ParseTermLocal := function(term_local)
        local parsed, idx_e1, idx_e2, idx_h1, idx_h2, idx_f1, idx_f2, root_local;
        if IsBoundGlobal("ParseTerm") then
            parsed := ValueGlobal("ParseTerm")(term_local);
            if IsRecord(parsed) and IsBound(parsed.kind) and IsBound(parsed.root) then
                return parsed;
            fi;
        fi;
        idx_e1 := PositionSublist(term_local, "E_{[");
        idx_e2 := PositionSublist(term_local, "]}");
        if idx_e1 <> fail and idx_e2 <> fail and idx_e2 > idx_e1 + 3 then
            root_local := term_local{[idx_e1+4..idx_e2-1]};
            return rec(kind := "E", root := root_local, coeff := 1);
        fi;
        idx_h1 := PositionSublist(term_local, "H_{");
        idx_h2 := PositionSublist(term_local, "}");
        if idx_h1 <> fail and idx_h2 <> fail and idx_h2 > idx_h1 + 2 then
            root_local := term_local{[idx_h1+3..idx_h2-1]};
            return rec(kind := "H", root := root_local, coeff := 1);
        fi;
        idx_f1 := PositionSublist(term_local, "F_{[");
        idx_f2 := PositionSublist(term_local, "]}");
        if idx_f1 <> fail and idx_f2 <> fail and idx_f2 > idx_f1 + 3 then
            root_local := term_local{[idx_f1+4..idx_f2-1]};
            return rec(kind := "F", root := root_local, coeff := 1);
        fi;
        return rec(kind := "", root := "", coeff := 0);
    end;
    SolveUniqueSystem := function(Ain, bin)
        local nloc, Ab, row, col, i, j, pivot, temp, factor, pivot_cols, rank, allzero, x, sumv;
        nloc := Length(Ain);
        if nloc = 0 or Length(bin) <> nloc then
            return fail;
        fi;
        Ab := [];
        for i in [1..nloc] do
            Ab[i] := ShallowCopy(Ain[i]);
            Add(Ab[i], bin[i]);
        od;
        pivot_cols := [];
        row := 1;
        for col in [1..nloc] do
            pivot := 0;
            for i in [row..nloc] do
                if Ab[i][col] <> 0 then
                    pivot := i;
                    break;
                fi;
            od;
            if pivot = 0 then
                continue;
            fi;
            if pivot <> row then
                temp := Ab[row];
                Ab[row] := Ab[pivot];
                Ab[pivot] := temp;
            fi;
            for i in [row+1..nloc] do
                factor := Ab[i][col] / Ab[row][col];
                for j in [col..nloc+1] do
                    Ab[i][j] := Ab[i][j] - factor * Ab[row][j];
                od;
            od;
            Add(pivot_cols, col);
            row := row + 1;
            if row > nloc then
                break;
            fi;
        od;
        rank := Length(pivot_cols);
        for i in [rank+1..nloc] do
            allzero := true;
            for j in [1..nloc] do
                if Ab[i][j] <> 0 then
                    allzero := false;
                    break;
                fi;
            od;
            if allzero and Ab[i][nloc+1] <> 0 then
                return fail;
            fi;
        od;
        if rank <> nloc then
            return fail;
        fi;
        x := List([1..nloc], i -> 0);
        for i in [rank, rank-1..1] do
            col := pivot_cols[i];
            sumv := Ab[i][nloc+1];
            for j in [col+1..nloc] do
                sumv := sumv - Ab[i][j] * x[j];
            od;
            x[col] := sumv / Ab[i][col];
        od;
        return x;
    end;
    if not IsRecord(triple) then
        return rec(ok := false, h := fail, e := fail, f := fail);
    fi;
    parser := fail;
    if IsBoundGlobal("ParseRootStringFlexible") then
        parser := ValueGlobal("ParseRootStringFlexible");
    elif IsBoundGlobal("ParseRootString") then
        parser := ValueGlobal("ParseRootString");
    fi;
    if parser = fail or not IsBoundGlobal("F4_SymMatrix") then
        return rec(ok := false, h := triple.h, e := triple.e, f := triple.f);
    fi;
    S := ValueGlobal("F4_SymMatrix")();
    e_terms := [];
    seen := [];
    if triple.h = "0" then
        return rec(ok := false, h := triple.h, e := triple.e, f := triple.f);
    fi;
    e_roots := [];
    for term in SplitByPlus(triple.e) do
        pair := ParseTermLocal(term);
        if pair.kind = "E" then
            if not pair.root in seen then
                Add(seen, pair.root);
                Add(e_roots, pair.root);
            fi;
        fi;
    od;
    if Length(e_roots) = 0 then
        for term in SplitByPlus(triple.h) do
            pair := ParseTermLocal(term);
            if pair.kind = "H" then
                if not pair.root in seen then
                    Add(seen, pair.root);
                    Add(e_roots, pair.root);
                fi;
            fi;
        od;
    fi;
    for term in e_roots do
        v := parser(term);
        if not IsList(v) or Length(v) <> 4 then
            return rec(ok := false, h := triple.h, e := triple.e, f := triple.f);
        fi;
        Add(e_terms, term);
    od;
    BuildSystemAndSolve := function(root_names)
        local roots_loc, nloc, Aloc, bloc, iloc, jloc, viloc, vjloc, numloc, denloc, cloc, okloc, lhsloc;
        roots_loc := [];
        for term in root_names do
            v := parser(term);
            if not IsList(v) or Length(v) <> 4 then
                return fail;
            fi;
            Add(roots_loc, v);
        od;
        nloc := Length(roots_loc);
        if nloc = 0 then
            return fail;
        fi;
        Aloc := [];
        bloc := [];
        for iloc in [1..nloc] do
            Aloc[iloc] := [];
            viloc := roots_loc[iloc];
            for jloc in [1..nloc] do
                vjloc := roots_loc[jloc];
                numloc := viloc[1]*(S[1][1]*vjloc[1]+S[1][2]*vjloc[2]+S[1][3]*vjloc[3]+S[1][4]*vjloc[4])
                       + viloc[2]*(S[2][1]*vjloc[1]+S[2][2]*vjloc[2]+S[2][3]*vjloc[3]+S[2][4]*vjloc[4])
                       + viloc[3]*(S[3][1]*vjloc[1]+S[3][2]*vjloc[2]+S[3][3]*vjloc[3]+S[3][4]*vjloc[4])
                       + viloc[4]*(S[4][1]*vjloc[1]+S[4][2]*vjloc[2]+S[4][3]*vjloc[3]+S[4][4]*vjloc[4]);
                denloc := vjloc[1]*(S[1][1]*vjloc[1]+S[1][2]*vjloc[2]+S[1][3]*vjloc[3]+S[1][4]*vjloc[4])
                       + vjloc[2]*(S[2][1]*vjloc[1]+S[2][2]*vjloc[2]+S[2][3]*vjloc[3]+S[2][4]*vjloc[4])
                       + vjloc[3]*(S[3][1]*vjloc[1]+S[3][2]*vjloc[2]+S[3][3]*vjloc[3]+S[3][4]*vjloc[4])
                       + vjloc[4]*(S[4][1]*vjloc[1]+S[4][2]*vjloc[2]+S[4][3]*vjloc[3]+S[4][4]*vjloc[4]);
                if denloc = 0 then
                    return fail;
                fi;
                Add(Aloc[iloc], 2 * numloc / denloc);
            od;
            Add(bloc, 2);
        od;
        cloc := SolveUniqueSystem(Aloc, bloc);
        if cloc = fail then
            return fail;
        fi;
        okloc := true;
        for iloc in [1..nloc] do
            lhsloc := 0;
            for jloc in [1..nloc] do
                lhsloc := lhsloc + Aloc[iloc][jloc] * cloc[jloc];
            od;
            if IsFloat(lhsloc) then
                if AbsFloat(lhsloc - 2) > 1.0e-6 then
                    okloc := false;
                    break;
                fi;
            else
                if lhsloc <> 2 then
                    okloc := false;
                    break;
                fi;
            fi;
        od;
        if not okloc then
            return fail;
        fi;
        return rec(c := cloc, A := Aloc, b := bloc, roots := roots_loc);
    end;
    ValidateAgainstERoots := function(alpha_names, gamma_names, coeffs)
        local ii, jj, alpha_v, gamma_v, lhs_v, num_v, den_v;
        if Length(gamma_names) <> Length(coeffs) then
            return false;
        fi;
        for ii in [1..Length(alpha_names)] do
            alpha_v := parser(alpha_names[ii]);
            if not IsList(alpha_v) or Length(alpha_v) <> 4 then
                return false;
            fi;
            lhs_v := 0;
            for jj in [1..Length(gamma_names)] do
                gamma_v := parser(gamma_names[jj]);
                if not IsList(gamma_v) or Length(gamma_v) <> 4 then
                    return false;
                fi;
                num_v := alpha_v[1]*(S[1][1]*gamma_v[1]+S[1][2]*gamma_v[2]+S[1][3]*gamma_v[3]+S[1][4]*gamma_v[4])
                      + alpha_v[2]*(S[2][1]*gamma_v[1]+S[2][2]*gamma_v[2]+S[2][3]*gamma_v[3]+S[2][4]*gamma_v[4])
                      + alpha_v[3]*(S[3][1]*gamma_v[1]+S[3][2]*gamma_v[2]+S[3][3]*gamma_v[3]+S[3][4]*gamma_v[4])
                      + alpha_v[4]*(S[4][1]*gamma_v[1]+S[4][2]*gamma_v[2]+S[4][3]*gamma_v[3]+S[4][4]*gamma_v[4]);
                den_v := gamma_v[1]*(S[1][1]*gamma_v[1]+S[1][2]*gamma_v[2]+S[1][3]*gamma_v[3]+S[1][4]*gamma_v[4])
                      + gamma_v[2]*(S[2][1]*gamma_v[1]+S[2][2]*gamma_v[2]+S[2][3]*gamma_v[3]+S[2][4]*gamma_v[4])
                      + gamma_v[3]*(S[3][1]*gamma_v[1]+S[3][2]*gamma_v[2]+S[3][3]*gamma_v[3]+S[3][4]*gamma_v[4])
                      + gamma_v[4]*(S[4][1]*gamma_v[1]+S[4][2]*gamma_v[2]+S[4][3]*gamma_v[3]+S[4][4]*gamma_v[4]);
                if den_v = 0 then
                    return false;
                fi;
                lhs_v := lhs_v + coeffs[jj] * (2 * num_v / den_v);
            od;
            if IsFloat(lhs_v) then
                if AbsFloat(lhs_v - 2) > 1.0e-6 then
                    return false;
                fi;
            else
                if lhs_v <> 2 then
                    return false;
                fi;
            fi;
        od;
        return true;
    end;
    n := Length(e_terms);
    if n = 0 then
        return rec(ok := false, h := triple.h, e := triple.e, f := triple.f);
    fi;
    c_try := BuildSystemAndSolve(e_terms);
    if c_try <> fail then
        c := c_try.c;
        A := c_try.A;
        b := c_try.b;
        roots := c_try.roots;
        validated_ok := ValidateAgainstERoots(e_roots, e_terms, c);
        ok := validated_ok;
    else
        ok := false;
    fi;
    if not ok then
        best_roots := fail;
        best_try := fail;
        chosen_roots := [];
        DFSSubsetSearch := function(start_idx, remaining)
            local idx, trial_roots, trial_try;
            if best_roots <> fail then
                return;
            fi;
            if remaining = 0 then
                trial_roots := Compacted(ShallowCopy(chosen_roots));
                trial_try := BuildSystemAndSolve(trial_roots);
                if trial_try <> fail and ValidateAgainstERoots(trial_roots, trial_roots, trial_try.c) then
                    best_roots := trial_roots;
                    best_try := trial_try;
                fi;
                return;
            fi;
            for idx in [start_idx..Length(e_terms) - remaining + 1] do
                Add(chosen_roots, e_terms[idx]);
                DFSSubsetSearch(idx + 1, remaining - 1);
                Unbind(chosen_roots[Length(chosen_roots)]);
                if best_roots <> fail then
                    return;
                fi;
            od;
        end;
        for subset_size in [Length(e_terms) - 1, Length(e_terms) - 2..1] do
            DFSSubsetSearch(1, subset_size);
            if best_roots <> fail then
                break;
            fi;
        od;
        if best_roots = fail or best_try = fail then
            return rec(ok := false, h := triple.h, e := triple.e, f := triple.f);
        fi;
        e_terms := best_roots;
        e_roots := ShallowCopy(best_roots);
        c := best_try.c;
        A := best_try.A;
        b := best_try.b;
        roots := best_try.roots;
        n := Length(e_terms);
    fi;
    h_expr := "";
    e_expr := "";
    f_expr := "";
    for i in [1..n] do
        coeff := c[i];
        if h_expr <> "" then h_expr := Concatenation(h_expr, " + "); fi;
        if coeff = 1 then
            h_expr := Concatenation(h_expr, "H_{", e_terms[i], "}");
        else
            h_expr := Concatenation(h_expr, String(coeff), "H_{", e_terms[i], "}");
        fi;
        if e_expr <> "" then e_expr := Concatenation(e_expr, " + "); fi;
        e_expr := Concatenation(e_expr, "E_{[", e_terms[i], "]}");
        if f_expr <> "" then f_expr := Concatenation(f_expr, " + "); fi;
        if coeff = 1 then
            f_expr := Concatenation(f_expr, "F_{[", e_terms[i], "]}");
        else
            f_expr := Concatenation(f_expr, String(coeff), "F_{[", e_terms[i], "]}");
        fi;
    od;
    if h_expr = "" then h_expr := "0"; fi;
    if e_expr = "" then e_expr := "0"; fi;
    if f_expr = "" then f_expr := "0"; fi;
    return rec(ok := true, h := h_expr, e := e_expr, f := f_expr);
end;;

RefitCompositeTripleF4_PostCheck := function(h_str, e_str, f_str, normalized_roots, effective_colors, efh_ok, pk_ok)
    local triple_refit, norm, term_f4, pair_f4, vec_f4, nc_h, nc_e, nc_f;
    if (not efh_ok or not pk_ok) and IsBoundGlobal("RefitCompositeTripleForF4") then
        triple_refit := ValueGlobal("RefitCompositeTripleForF4")(rec(h := h_str, e := e_str, f := f_str));
        if IsRecord(triple_refit) and IsBound(triple_refit.ok) and triple_refit.ok = true then
            h_str := triple_refit.h;
            e_str := triple_refit.e;
            f_str := triple_refit.f;
            Print("        [重整] 复合项失败后强制重解: PASS\n");
            if IsBoundGlobal("VerifyEFH_NormalConditions") then
                efh_ok := ValueGlobal("VerifyEFH_NormalConditions")(h_str, e_str, f_str);
            fi;
            if IsBoundGlobal("VerifyPK_Subsystem_F4") then
                pk_ok := ValueGlobal("VerifyPK_Subsystem_F4")(h_str, e_str, f_str, normalized_roots, effective_colors);
            fi;
            if (not efh_ok or not pk_ok) and IsBoundGlobal("NormalizeSL2Triple_Unit") then
                norm := ValueGlobal("NormalizeSL2Triple_Unit")(h_str, e_str, f_str);
                if IsRecord(norm) and IsBound(norm.ok) and norm.ok = true then
                    triple_refit := ValueGlobal("RefitCompositeTripleForF4")(rec(h := norm.h, e := norm.e, f := norm.f));
                    if IsRecord(triple_refit) and IsBound(triple_refit.ok) and triple_refit.ok = true then
                        h_str := triple_refit.h;
                        e_str := triple_refit.e;
                        f_str := triple_refit.f;
                        Print("        [重整] 复合项失败后二次重解: PASS\n");
                        if IsBoundGlobal("VerifyEFH_NormalConditions") then
                            efh_ok := ValueGlobal("VerifyEFH_NormalConditions")(h_str, e_str, f_str);
                        fi;
                        if IsBoundGlobal("VerifyPK_Subsystem_F4") then
                            pk_ok := ValueGlobal("VerifyPK_Subsystem_F4")(h_str, e_str, f_str, normalized_roots, effective_colors);
                        fi;
                    else
                        Print("        [重整] 复合项失败后二次重解: FAIL\n");
                    fi;
                fi;
            fi;
            if (not efh_ok or not pk_ok) and IsBoundGlobal("ParseRootStringFlexible") and IsBoundGlobal("IsNoncompactF4Vec") then
                nc_h := "";
                nc_e := "";
                nc_f := "";
                for term_f4 in SplitString(e_str, " + ") do
                    pair_f4 := ValueGlobal("ParseTerm")(term_f4);
                    if IsRecord(pair_f4) and IsBound(pair_f4.kind) and IsBound(pair_f4.root) and pair_f4.kind = "E" then
                        vec_f4 := ValueGlobal("ParseRootStringFlexible")(pair_f4.root);
                        if IsList(vec_f4) and Length(vec_f4) = 4 and ValueGlobal("IsNoncompactF4Vec")(vec_f4) then
                            if nc_h <> "" then nc_h := Concatenation(nc_h, " + "); fi;
                            nc_h := Concatenation(nc_h, "H_{", pair_f4.root, "}");
                            if nc_e <> "" then nc_e := Concatenation(nc_e, " + "); fi;
                            nc_e := Concatenation(nc_e, "E_{[", pair_f4.root, "]}");
                            if nc_f <> "" then nc_f := Concatenation(nc_f, " + "); fi;
                            nc_f := Concatenation(nc_f, "F_{[", pair_f4.root, "]}");
                        fi;
                    fi;
                od;
                if nc_e <> "" then
                    triple_refit := ValueGlobal("RefitCompositeTripleForF4")(rec(h := nc_h, e := nc_e, f := nc_f));
                    if IsRecord(triple_refit) and IsBound(triple_refit.ok) and triple_refit.ok = true then
                        h_str := triple_refit.h;
                        e_str := triple_refit.e;
                        f_str := triple_refit.f;
                        Print("        [重整] 复合项非紧过滤重解: PASS\n");
                        if IsBoundGlobal("VerifyEFH_NormalConditions") then
                            efh_ok := ValueGlobal("VerifyEFH_NormalConditions")(h_str, e_str, f_str);
                        fi;
                        if IsBoundGlobal("VerifyPK_Subsystem_F4") then
                            pk_ok := ValueGlobal("VerifyPK_Subsystem_F4")(h_str, e_str, f_str, normalized_roots, effective_colors);
                        fi;
                    else
                        Print("        [重整] 复合项非紧过滤重解: FAIL\n");
                    fi;
                fi;
            fi;
        else
            Print("        [重整] 复合项失败后强制重解: FAIL\n");
        fi;
    fi;
    return rec(h := h_str, e := e_str, f := f_str, efh_ok := efh_ok, pk_ok := pk_ok);
end;;

NormalizeThetaCandidateTripleForF4 := function(triple, normalized_roots, effective_colors)
    local out, efh_ok, pk_ok, triple_refit;
    if not IsRecord(triple) then
        return rec(triple := triple, efh_ok := false, pk_ok := false);
    fi;
    out := ShallowCopy(triple);
    efh_ok := false;
    if IsBoundGlobal("VerifyEFH_NormalConditions") then
        efh_ok := ValueGlobal("VerifyEFH_NormalConditions")(out.h, out.e, out.f);
    fi;
    pk_ok := true;
    if IsBoundGlobal("VerifyPK_Subsystem_F4") then
        pk_ok := ValueGlobal("VerifyPK_Subsystem_F4")(out.h, out.e, out.f, normalized_roots, effective_colors);
    fi;
    if (not efh_ok or not pk_ok) and IsBoundGlobal("RefitCompositeTripleForF4") then
        triple_refit := ValueGlobal("RefitCompositeTripleForF4")(rec(h := out.h, e := out.e, f := out.f));
        if IsRecord(triple_refit) and IsBound(triple_refit.ok) and triple_refit.ok = true then
            out.h := triple_refit.h;
            out.e := triple_refit.e;
            out.f := triple_refit.f;
            if IsBoundGlobal("VerifyEFH_NormalConditions") then
                efh_ok := ValueGlobal("VerifyEFH_NormalConditions")(out.h, out.e, out.f);
            fi;
            pk_ok := true;
            if IsBoundGlobal("VerifyPK_Subsystem_F4") then
                pk_ok := ValueGlobal("VerifyPK_Subsystem_F4")(out.h, out.e, out.f, normalized_roots, effective_colors);
            fi;
        fi;
    fi;
    return rec(triple := out, efh_ok := efh_ok, pk_ok := pk_ok);
end;;

DispatchKOrbitOutput := function(ambient, type, h_str, e_str, f_str)
    if ambient = "G2" and IsBoundGlobal("PrintKOrbitInfo") then
        GlobalCurrentSubsystemType := type;
        ValueGlobal("PrintKOrbitInfo")(h_str);
        return;
    fi;
    if IsBoundGlobal("PrintKOrbitInfo_F4") then
        if IsBoundGlobal("GlobalColorsListForDiagram") then
            if Length(ValueGlobal("GlobalColorsListForDiagram")) = 4 then
                ValueGlobal("PrintKOrbitInfo_F4")(rec(h := h_str, e := e_str, f := f_str), Concatenation(ValueGlobal("GlobalColorsListForDiagram"), [0]));
            elif Length(ValueGlobal("GlobalColorsListForDiagram")) >= 5 then
                ValueGlobal("PrintKOrbitInfo_F4")(rec(h := h_str, e := e_str, f := f_str), ValueGlobal("GlobalColorsListForDiagram"));
            else
                ValueGlobal("PrintKOrbitInfo_F4")(rec(h := h_str, e := e_str, f := f_str), [0,0,0,0,0]);
            fi;
        else
            ValueGlobal("PrintKOrbitInfo_F4")(rec(h := h_str, e := e_str, f := f_str), [0,0,0,0,0]);
        fi;
        return;
    fi;
    if IsBoundGlobal("PrintKOrbitInfo") then
        ValueGlobal("PrintKOrbitInfo")(h_str);
    else
        ValueGlobal("PrintGenericHSummary")(h_str);
    fi;
end;;

NormalizeExprForSplitTermsF4 := function(expr)
    local s;
    s := Filtered(expr, c -> c <> '\n' and c <> '\r' and c <> '\t' and c <> '\\');
    while PositionSublist(s, "  ") <> fail do
        s := ReplacedString(s, "  ", " ");
    od;
    return s;
end;;

SplitExprTermsByPlusF4 := function(expr)
    local cleaned, terms, rest, pos;
    cleaned := NormalizeExprForSplitTermsF4(expr);
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

F4TripleFilterInnerProduct := function(a, b)
    return a[1] * (4 * b[1] - 2 * b[2]) +
           a[2] * (-2 * b[1] + 4 * b[2] - 2 * b[3]) +
           a[3] * (-2 * b[2] + 2 * b[3] - b[4]) +
           a[4] * (-b[3] + 2 * b[4]);
end;;

F4TripleFilterReflectByBeta := function(v, beta)
    local beta_len2, coeff, i, out;
    beta_len2 := F4TripleFilterInnerProduct(beta, beta);
    coeff := 2 * F4TripleFilterInnerProduct(v, beta) / beta_len2;
    out := [];
    for i in [1..4] do
        Add(out, v[i] - coeff * beta[i]);
    od;
    return out;
end;;

F4TripleFilterReflectBySimpleIndex := function(v, idx)
    local beta;
    if idx = 0 then
        beta := [-2, -3, -4, -2];
    else
        beta := [0, 0, 0, 0];
        beta[idx] := 1;
    fi;
    return F4TripleFilterReflectByBeta(v, beta);
end;;

F4TripleFilterReflectByRoot := function(v, beta)
    return F4TripleFilterReflectByBeta(v, beta);
end;;

GenerateDeltaJRootSetForFilterF4 := function(simple_roots)
    local roots, queue, seen, beta, root, next_root, key_of, format_root;
    if not IsBoundGlobal("FormatRootF4") then
        return [];
    fi;
    format_root := ValueGlobal("FormatRootF4");
    key_of := v -> format_root(v);
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
            next_root := F4TripleFilterReflectByRoot(root, beta);
            if Position(seen, key_of(next_root)) = fail then
                Add(roots, next_root);
                Add(queue, next_root);
                Add(seen, key_of(next_root));
            fi;
        od;
    od;
    return roots;
end;;

ApplyWKWordOnRootF4ForFilter := function(v, word)
    local out, idx;
    out := [v[1], v[2], v[3], v[4]];
    for idx in word do
        out := F4TripleFilterReflectBySimpleIndex(out, idx);
    od;
    return out;
end;;

BuildConjugacyKeyForExprF4 := function(expr, expected_kind, word, parser)
    local terms, t, data, roots, v, wv, parse_term, format_root;
    if not IsBoundGlobal("ParseTerm") or not IsBoundGlobal("FormatRootF4") then
        return fail;
    fi;
    parse_term := ValueGlobal("ParseTerm");
    format_root := ValueGlobal("FormatRootF4");
    terms := SplitExprTermsByPlusF4(expr);
    roots := [];
    for t in terms do
        data := parse_term(t);
        if data.kind = expected_kind then
            v := parser(data.root);
            if v = fail or not IsList(v) or Length(v) <> 4 then
                return fail;
            fi;
            wv := ApplyWKWordOnRootF4ForFilter(v, word);
            Add(roots, Concatenation(String(data.coeff), "@", format_root(wv)));
        fi;
    od;
    return SortedList(roots);
end;;

IsWKConjugateFixingDeltaJF4 := function(t1, t2, subsystem_roots)
    local parser, roots_vecs, r, delta_roots, delta_key, generators, queue, seen_states, state, gen, next_state, key_state, key2_h, key2_e, key2_f, key1_h, key1_e, key1_f, key_of_roots, basis0, basis_key, ApplyGenOnBasis, stabilizer_words, moved, format_root;
    parser := fail;
    if IsBoundGlobal("ParseRootStringFlexible") then
        parser := ValueGlobal("ParseRootStringFlexible");
    elif IsBoundGlobal("ParseRootString") then
        parser := ValueGlobal("ParseRootString");
    fi;
    if parser = fail then
        return false;
    fi;
    if not IsBoundGlobal("FormatRootF4") then
        return false;
    fi;
    format_root := ValueGlobal("FormatRootF4");
    roots_vecs := [];
    for r in subsystem_roots do
        if IsList(r) and Length(r) = 4 then
            Add(roots_vecs, [r[1], r[2], r[3], r[4]]);
        else
            Add(roots_vecs, parser(r));
        fi;
    od;
    if Length(roots_vecs) = 0 then
        return false;
    fi;
    delta_roots := GenerateDeltaJRootSetForFilterF4(roots_vecs);
    if Length(delta_roots) = 0 then
        return false;
    fi;
    key_of_roots := function(vecs)
        return SortedList(List(vecs, v -> format_root(v)));
    end;
    delta_key := key_of_roots(delta_roots);
    generators := [0, 2, 3, 4];
    basis0 := [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
    basis_key := function(images)
        return Concatenation(format_root(images[1]), "|", format_root(images[2]), "|", format_root(images[3]), "|", format_root(images[4]));
    end;
    ApplyGenOnBasis := function(images, idx)
        return List(images, v -> F4TripleFilterReflectBySimpleIndex(v, idx));
    end;
    state := rec(word := [], images := basis0);
    queue := [state];
    seen_states := [basis_key(basis0)];
    stabilizer_words := [];
    while Length(queue) > 0 do
        state := queue[1];
        Remove(queue, 1);
        moved := key_of_roots(List(delta_roots, x -> ApplyWKWordOnRootF4ForFilter(x, state.word)));
        if moved = delta_key then
            Add(stabilizer_words, state.word);
        fi;
        for gen in generators do
            next_state := rec(
                word := Concatenation(state.word, [gen]),
                images := ApplyGenOnBasis(state.images, gen)
            );
            key_state := basis_key(next_state.images);
            if Position(seen_states, key_state) = fail then
                Add(seen_states, key_state);
                Add(queue, next_state);
            fi;
        od;
    od;
    if Length(stabilizer_words) = 0 then
        return false;
    fi;
    key2_h := BuildConjugacyKeyForExprF4(t2.h, "H", [], parser);
    key2_e := BuildConjugacyKeyForExprF4(t2.e, "E", [], parser);
    key2_f := BuildConjugacyKeyForExprF4(t2.f, "F", [], parser);
    if key2_h = fail or key2_e = fail or key2_f = fail then
        return false;
    fi;
    for state in stabilizer_words do
        key1_h := BuildConjugacyKeyForExprF4(t1.h, "H", state, parser);
        key1_e := BuildConjugacyKeyForExprF4(t1.e, "E", state, parser);
        key1_f := BuildConjugacyKeyForExprF4(t1.f, "F", state, parser);
        if key1_e <> fail and key1_f <> fail and key1_e = key2_e and key1_f = key2_f then
            return true;
        fi;
    od;
    return false;
end;;

FilterTriplesByWKConjugacyFixingDeltaJF4 := function(triples, subsystem_roots)
    local kept, i, j, is_dup;
    kept := [];
    for i in [1..Length(triples)] do
        is_dup := false;
        for j in [1..Length(kept)] do
            if IsWKConjugateFixingDeltaJF4(triples[i], kept[j], subsystem_roots) then
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

FilterCompositeTripleFamiliesByKOrbitLabelF4 := function(triples)
    local calc_fn, path_len, i, item, path, calc, label, KeyOfPathSlice, JoinEntries, current_level, current_nodes, kept_nodes, seen_signatures, node, signature, parent_nodes, parent_path, parent_key, pos, child_sigs, kept_items;
    if not IsList(triples) or Length(triples) = 0 then
        return triples;
    fi;
    if not IsBoundGlobal("F4ComputeKCDominantFromTriple") then
        return triples;
    fi;
    calc_fn := ValueGlobal("F4ComputeKCDominantFromTriple");
    path_len := fail;
    KeyOfPathSlice := function(one_path, first_idx, last_idx)
        local out, idx;
        if first_idx > last_idx then
            return "";
        fi;
        out := String(one_path[first_idx]);
        for idx in [first_idx + 1..last_idx] do
            out := Concatenation(out, ",", String(one_path[idx]));
        od;
        return out;
    end;
    JoinEntries := function(entries)
        local out, idx;
        if Length(entries) = 0 then
            return "";
        fi;
        out := entries[1];
        for idx in [2..Length(entries)] do
            out := Concatenation(out, ";", entries[idx]);
        od;
        return out;
    end;
    current_nodes := [];
    for i in [1..Length(triples)] do
        item := triples[i];
        if not IsRecord(item) or not IsBound(item.component_orbit_choice) or not IsList(item.component_orbit_choice) or Length(item.component_orbit_choice) < 2 then
            return triples;
        fi;
        path := item.component_orbit_choice;
        if path_len = fail then
            path_len := Length(path);
        elif Length(path) <> path_len then
            return triples;
        fi;
        calc := calc_fn(rec(h := item.h, e := item.e, f := item.f));
        if not IsRecord(calc) or not IsBound(calc.ok) or not calc.ok or not IsBound(calc.resolved_label) then
            return triples;
        fi;
        label := calc.resolved_label;
        Add(current_nodes, rec(
            path := ShallowCopy(path),
            parent_key := KeyOfPathSlice(path, 1, path_len - 1),
            sig := label,
            items := [item]
        ));
    od;
    current_level := path_len;
    while current_level >= 1 do
        seen_signatures := [];
        kept_nodes := [];
        for node in current_nodes do
            signature := Concatenation(node.parent_key, "||", node.sig);
            if Position(seen_signatures, signature) = fail then
                Add(seen_signatures, signature);
                Add(kept_nodes, node);
            fi;
        od;
        if current_level = 1 then
            kept_items := [];
            for node in kept_nodes do
                Append(kept_items, node.items);
            od;
            return kept_items;
        fi;
        parent_nodes := [];
        for node in kept_nodes do
            parent_path := node.path{[1..current_level - 1]};
            parent_key := KeyOfPathSlice(parent_path, 1, current_level - 2);
            pos := PositionProperty(parent_nodes, g -> g.path = parent_path);
            if pos = fail then
                Add(parent_nodes, rec(
                    path := parent_path,
                    parent_key := parent_key,
                    child_sigs := [],
                    items := []
                ));
                pos := Length(parent_nodes);
            fi;
            child_sigs := parent_nodes[pos].child_sigs;
            Add(child_sigs, node.sig);
            parent_nodes[pos].child_sigs := child_sigs;
            Append(parent_nodes[pos].items, node.items);
        od;
        for i in [1..Length(parent_nodes)] do
            child_sigs := SortedList(parent_nodes[i].child_sigs);
            parent_nodes[i].sig := JoinEntries(child_sigs);
            Unbind(parent_nodes[i].child_sigs);
        od;
        current_nodes := parent_nodes;
        current_level := current_level - 1;
    od;
    return triples;
end;;

ResolveModule4RealFormInfo := function(type, rank)
    local info;
    if IsBoundGlobal("GlobalModule2CurrentRealFormInfo") then
        info := ValueGlobal("GlobalModule2CurrentRealFormInfo");
        if IsRecord(info) and IsBound(info.type) and IsBound(info.params) then
            return info;
        fi;
    fi;
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
    local rf_type, rank_a, rank_c, auts, p, q, theta, g, dims, target_k_dim, target_p_dim, finite_order_inner_automorphisms,
          total_dim, series_type, series_rank;
    if real_form_info = fail or not IsRecord(real_form_info) then
        return fail;
    fi;
    if not (IsBoundGlobal("LoadPackage") and LoadPackage("sla") = true) then
        return fail;
    fi;
    if not IsBoundGlobal("FiniteOrderInnerAutomorphisms") then
        return fail;
    fi;
    finite_order_inner_automorphisms := ValueGlobal("FiniteOrderInnerAutomorphisms");
    if not IsBound(real_form_info.type) then
        return fail;
    fi;
    rf_type := real_form_info.type;
    if rf_type = "su_pq" and IsBound(real_form_info.params) and Length(real_form_info.params) >= 2 then
        p := real_form_info.params[1];
        q := real_form_info.params[2];
        if IsInt(p) and IsInt(q) and p > 0 and q > 0 then
            rank_a := p + q - 1;
            if IsInt(rank_a) and rank_a >= 1 then
                auts := finite_order_inner_automorphisms("A", rank_a, 2);
                if IsList(auts) then
                    target_k_dim := p * p + q * q - 1;
                    target_p_dim := 2 * p * q;
                    for theta in auts do
                        g := Grading(theta);
                        if IsList(g) and Length(g) >= 2 and IsList(g[1]) and IsList(g[2]) then
                            dims := [Length(g[1]), Length(g[2])];
                            if dims[1] = target_k_dim and dims[2] = target_p_dim then
                                return theta;
                            fi;
                        fi;
                    od;
                    if Length(auts) > 0 then
                        return auts[1];
                    fi;
                fi;
            fi;
        fi;
    fi;
    if rf_type = "sl_R" and IsBound(real_form_info.params) and Length(real_form_info.params) >= 1 then
        rank_a := real_form_info.params[1];
        if IsInt(rank_a) and rank_a >= 1 then
            auts := finite_order_inner_automorphisms("A", rank_a, 2);
            if IsList(auts) and Length(auts) > 0 then
                return auts[1];
            fi;
        fi;
    fi;
    if rf_type = "sp_R" and IsBound(real_form_info.params) and Length(real_form_info.params) >= 1 then
        rank_c := real_form_info.params[1];
        if IsInt(rank_c) and rank_c >= 1 then
            auts := finite_order_inner_automorphisms("C", rank_c, 2);
            if IsList(auts) then
                target_k_dim := rank_c * rank_c;
                target_p_dim := rank_c * (rank_c + 1);
                for theta in auts do
                    g := Grading(theta);
                    if IsList(g) and Length(g) >= 2 and IsList(g[1]) and IsList(g[2]) then
                        dims := [Length(g[1]), Length(g[2])];
                        if dims[1] = target_k_dim and dims[2] = target_p_dim then
                            return theta;
                        fi;
                    fi;
                od;
                if Length(auts) > 0 then
                    return auts[1];
                fi;
            fi;
        fi;
    fi;
    if rf_type = "sp_pq" and IsBound(real_form_info.params) and Length(real_form_info.params) >= 2 then
        p := real_form_info.params[1];
        q := real_form_info.params[2];
        if IsInt(p) and IsInt(q) and p >= 0 and q >= 0 then
            rank_c := p + q;
            if rank_c >= 1 then
                auts := finite_order_inner_automorphisms("C", rank_c, 2);
                if IsList(auts) then
                    target_k_dim := p * (2 * p + 1) + q * (2 * q + 1);
                    target_p_dim := 4 * p * q;
                    for theta in auts do
                        g := Grading(theta);
                        if IsList(g) and Length(g) >= 2 and IsList(g[1]) and IsList(g[2]) then
                            dims := [Length(g[1]), Length(g[2])];
                            if dims[1] = target_k_dim and dims[2] = target_p_dim then
                                return theta;
                            fi;
                        fi;
                    od;
                fi;
            fi;
        fi;
    fi;
    if rf_type = "so_pq" and IsBound(real_form_info.params) and Length(real_form_info.params) >= 2 then
        p := real_form_info.params[1];
        q := real_form_info.params[2];
        if IsInt(p) and IsInt(q) and p >= 0 and q >= 0 then
            total_dim := p + q;
            if total_dim >= 3 then
                if (total_dim mod 2) = 1 then
                    series_type := "B";
                    series_rank := (total_dim - 1) / 2;
                else
                    series_type := "D";
                    series_rank := total_dim / 2;
                fi;
                auts := finite_order_inner_automorphisms(series_type, series_rank, 2);
                if IsList(auts) then
                    target_k_dim := p * (p - 1) / 2 + q * (q - 1) / 2;
                    target_p_dim := p * q;
                    for theta in auts do
                        g := Grading(theta);
                        if IsList(g) and Length(g) >= 2 and IsList(g[1]) and IsList(g[2]) then
                            dims := [Length(g[1]), Length(g[2])];
                            if dims[1] = target_k_dim and dims[2] = target_p_dim then
                                return theta;
                            fi;
                        fi;
                    od;
                fi;
            fi;
        fi;
    fi;
    return fail;
end;;

BuildThetaRepresentationInputFromRealForm := function(real_form_info)
    local theta;
    if real_form_info = fail or not IsRecord(real_form_info) then
        return fail;
    fi;
    if not (IsBoundGlobal("LoadPackage") and LoadPackage("sla") = true) then
        return fail;
    fi;
    if not IsBound(real_form_info.type) then
        return fail;
    fi;
    theta := BuildThetaAutomorphismFromRealForm(real_form_info);
    if theta <> fail then
        return rec(
            ok := true,
            mode := "automorphism",
            real_form := real_form_info,
            theta := theta
        );
    fi;
    return fail;
end;;

ComputeThetaOrbitTriplesForModule4 := function(real_form_info)
    local theta_input, theta, orbs, triples, i, t, cdim, noticed, selected, kac, zdim, L, g, l0_basis, lk, zlk, e_subalg, c_full, c, zc, triple_subalg, c_triple_full, c_triple, triple_cdim, zc_dim, noticed_triples, grading, nilpotent_orbits_of_theta_representation, has_kac_diagram, kac_diagram;
    if real_form_info = fail then
        return rec(ok := false, reason := "无实形式上下文");
    fi;
    theta_input := BuildThetaRepresentationInputFromRealForm(real_form_info);
    if theta_input = fail then
        return rec(ok := false, reason := "无法构造 theta/分次 数据");
    fi;
    if IsBoundGlobal("SetInfoLevel") and IsBoundGlobal("InfoSLA") then
        ValueGlobal("SetInfoLevel")(ValueGlobal("InfoSLA"), 1);
    fi;
    if not IsBound(theta_input.theta) or theta_input.theta = fail then
        return rec(ok := false, reason := "无法构造 theta 自同构");
    fi;
    if not IsBoundGlobal("NilpotentOrbitsOfThetaRepresentation") then
        return rec(ok := false, reason := "未加载 NilpotentOrbitsOfThetaRepresentation");
    fi;
    theta := theta_input.theta;
    grading := fail;
    nilpotent_orbits_of_theta_representation := ValueGlobal("NilpotentOrbitsOfThetaRepresentation");
    orbs := nilpotent_orbits_of_theta_representation(theta);
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
    if theta <> fail and IsBoundGlobal("HasKacDiagram") and IsBoundGlobal("KacDiagram") then
        has_kac_diagram := ValueGlobal("HasKacDiagram");
        kac_diagram := ValueGlobal("KacDiagram");
        if has_kac_diagram(theta) then
            kac := kac_diagram(theta).weights;
        else
            kac := fail;
        fi;
    else
        kac := fail;
    fi;
    return rec(
        ok := true,
        real_form := real_form_info,
        theta_input_mode := "automorphism",
        theta := theta,
        grading := grading,
        kac_weights := kac,
        center_dim_lk := zdim,
        triples := triples,
        noticed_triples := noticed_triples,
        selected := selected
    );
end;;

PrintThetaOrbitStep2 := function(theta_data)
    local t;
    if theta_data = fail or not IsRecord(theta_data) or not IsBound(theta_data.ok) or not theta_data.ok then
        if theta_data <> fail and IsRecord(theta_data) and IsBound(theta_data.reason) then
            Print("        [M4-STEP2] Theta-orbit 路径跳过: ", theta_data.reason, "\n");
        else
            Print("        [M4-STEP2] Theta-orbit 路径跳过\n");
        fi;
        return;
    fi;
    Print("        [M4-STEP2] NilpotentOrbitsOfThetaRepresentation\n");
    Print("        [M4-STEP2] RealForm = ", theta_data.real_form.type, ", Params = ", theta_data.real_form.params, "\n");
    if IsBound(theta_data.theta_input_mode) then
        Print("        [M4-STEP2] 输入模式 = ", theta_data.theta_input_mode, "\n");
    fi;
    if IsBound(theta_data.grading) and theta_data.grading <> fail then
        Print("        [M4-STEP2] grading d = ", theta_data.grading, "\n");
    fi;
    if IsBound(theta_data.kac_weights) and theta_data.kac_weights <> fail then
        Print("        [M4-STEP2] theta Kac weights = ", theta_data.kac_weights, "\n");
    fi;
    if IsBound(theta_data.center_dim_lk) and theta_data.center_dim_lk <> fail then
        Print("        [M4-STEP2] dim(z(l)∩k) = ", theta_data.center_dim_lk, "\n");
    fi;
    Print("        [M4-STEP2] theta 轨道总数 = ", Length(theta_data.triples), "\n");
    Print("        [M4-STEP2] noticed 轨道数 = ", Length(theta_data.noticed_triples), "\n");
    for t in theta_data.noticed_triples do
        Print("        [M4-STEP2] orbit#", t.index, ": x=", t.x, ", h=", t.h, ", y=", t.y, "\n");
        Print("        [M4-STEP2] orbit#", t.index, ": dim(C_e)= ", t.centralizer_dim_lk, ", dim(C_e∩Z)= ", t.center_intersection_dim, ", dim((l∩k)^[x,e,f])= ", t.triple_centralizer_dim_lk, ", dim(z(l)∩k)= ", t.center_dim_lk, ", noticed=", t.noticed, "\n");
    od;
    if Length(theta_data.noticed_triples) = 0 then
        Print("        [M4-STEP2] 未发现 noticed 轨道\n");
    fi;
    if theta_data.selected <> fail then
        Print("        [M4-STEP2] 选中用于根式输出的 triple: orbit#", theta_data.selected.index, " (noticed=", theta_data.selected.noticed, ")\n");
    fi;
end;;

NegateRootStringF4 := function(root_str)
    local parser, vec;
    parser := fail;
    if IsBoundGlobal("ParseRootStringFlexible") then
        parser := ValueGlobal("ParseRootStringFlexible");
    elif IsBoundGlobal("ParseRootString") then
        parser := ValueGlobal("ParseRootString");
    fi;
    if parser = fail then
        return root_str;
    fi;
    vec := parser(root_str);
    if not IsList(vec) then
        return root_str;
    fi;
    vec := List(vec, x -> -x);
    if IsBoundGlobal("FormatRootF4") then
        return ValueGlobal("FormatRootF4")(vec);
    fi;
    return FormatRoot(vec);
end;;

BuildUnitRootTripleFromList := function(root_names)
    local h_expr, e_expr, f_expr, idx;
    h_expr := "";
    e_expr := "";
    f_expr := "";
    for idx in [1..Length(root_names)] do
        if h_expr <> "" then
            h_expr := Concatenation(h_expr, " + ");
            e_expr := Concatenation(e_expr, " + ");
            f_expr := Concatenation(f_expr, " + ");
        fi;
        h_expr := Concatenation(h_expr, "H_{", root_names[idx], "}");
        e_expr := Concatenation(e_expr, "E_{[", root_names[idx], "]}");
        f_expr := Concatenation(f_expr, "F_{[", root_names[idx], "]}");
    od;
    if h_expr = "" then h_expr := "0"; fi;
    if e_expr = "" then e_expr := "0"; fi;
    if f_expr = "" then f_expr := "0"; fi;
    return rec(h := h_expr, e := e_expr, f := f_expr);
end;;

BuildDirectThetaRootListForF4 := function(theta_data, subsystem_roots, colors_opt)
    local selected, real_form, parser, black_positions, black_idx, root_vecs, root_name, vec, left_vec, mid_vec, right_vec, chosen_roots, negated;
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
    if not (real_form.type = "su_pq" and IsBound(real_form.params) and Length(real_form.params) >= 2 and real_form.params[1] = 2 and real_form.params[2] = 2) then
        return fail;
    fi;
    if not (IsList(subsystem_roots) and Length(subsystem_roots) = 3 and IsList(colors_opt) and Length(colors_opt) = 3) then
        return fail;
    fi;
    if IsBoundGlobal("ParseRootStringFlexible") then
        parser := ValueGlobal("ParseRootStringFlexible");
    elif IsBoundGlobal("ParseRootString") then
        parser := ValueGlobal("ParseRootString");
    else
        return fail;
    fi;
    black_positions := Filtered([1..Length(colors_opt)], i -> colors_opt[i] = 1);
    if Length(black_positions) <> 1 then
        return fail;
    fi;
    black_idx := black_positions[1];
    if black_idx <= 1 or black_idx >= Length(subsystem_roots) then
        return fail;
    fi;
    root_vecs := [];
    for root_name in subsystem_roots do
        vec := parser(root_name);
        if not (IsList(vec) and Length(vec) = 4) then
            return fail;
        fi;
        Add(root_vecs, vec);
    od;
    left_vec := [
        root_vecs[black_idx - 1][1] + root_vecs[black_idx][1],
        root_vecs[black_idx - 1][2] + root_vecs[black_idx][2],
        root_vecs[black_idx - 1][3] + root_vecs[black_idx][3],
        root_vecs[black_idx - 1][4] + root_vecs[black_idx][4]
    ];
    mid_vec := [
        -root_vecs[black_idx][1],
        -root_vecs[black_idx][2],
        -root_vecs[black_idx][3],
        -root_vecs[black_idx][4]
    ];
    right_vec := [
        root_vecs[black_idx][1] + root_vecs[black_idx + 1][1],
        root_vecs[black_idx][2] + root_vecs[black_idx + 1][2],
        root_vecs[black_idx][3] + root_vecs[black_idx + 1][3],
        root_vecs[black_idx][4] + root_vecs[black_idx + 1][4]
    ];
    if IsBoundGlobal("IsF4RootVec") then
        if not ValueGlobal("IsF4RootVec")(left_vec) or not ValueGlobal("IsF4RootVec")(mid_vec) or not ValueGlobal("IsF4RootVec")(right_vec) then
            return fail;
        fi;
    fi;
    chosen_roots := [
        ValueGlobal("FormatRootF4")(left_vec),
        ValueGlobal("FormatRootF4")(mid_vec),
        ValueGlobal("FormatRootF4")(right_vec)
    ];
    negated := false;
    if IsBound(selected.index) and ((selected.index mod 2) = 0) then
        chosen_roots := List(chosen_roots, r -> NegateRootStringF4(r));
        negated := true;
    fi;
    return rec(roots := chosen_roots, source_x := String(selected.x), source_mode := "theta_direct", negated := negated);
end;;

FindBasisSlotForLieVectorF4 := function(BL, vec)
    local coeffs, nz;
    coeffs := Coefficients(BL, vec);
    nz := Filtered([1..Length(coeffs)], i -> coeffs[i] <> 0);
    if Length(nz) <> 1 then
        return fail;
    fi;
    return rec(index := nz[1], scale := coeffs[nz[1]]);
end;;

GetScaledCoefficientFromSlotF4 := function(coeffs, slot)
    if slot = fail or not IsRecord(slot) or not IsBound(slot.index) or not IsBound(slot.scale) then
        return fail;
    fi;
    if slot.scale = 0 then
        return fail;
    fi;
    return coeffs[slot.index] / slot.scale;
end;;

FormatExprFromMapOrZeroF4 := function(map, prefix)
    if not IsList(map) or Length(map) = 0 then
        return "0";
    fi;
    if not IsBoundGlobal("FormatExprFromMap") then
        return "0";
    fi;
    return ValueGlobal("FormatExprFromMap")(map, prefix);
end;;

BuildThetaSourceRootDataF4 := function(theta_data)
    local L, R, S, Sinv, pos_roots, pos_vecs, neg_vecs, cartan_vecs, BL, pos_slots, neg_slots, cartan_slots, i;
    if theta_data = fail or not IsRecord(theta_data) or not IsBound(theta_data.theta) or theta_data.theta = fail then
        return fail;
    fi;
    if not IsBoundGlobal("RootSystem") or not IsBoundGlobal("SimpleRootsAsWeights")
       or not IsBoundGlobal("PositiveRootsAsWeights") or not IsBoundGlobal("PositiveRootVectors")
       or not IsBoundGlobal("NegativeRootVectors") or not IsBoundGlobal("CanonicalGenerators") then
        return fail;
    fi;
    L := Source(theta_data.theta);
    R := RootSystem(L);
    S := ValueGlobal("SimpleRootsAsWeights")(R);
    Sinv := S^-1;
    pos_roots := List(PositiveRootsAsWeights(R), w -> List(w * Sinv, x -> Int(x)));
    pos_vecs := PositiveRootVectors(R);
    neg_vecs := NegativeRootVectors(R);
    cartan_vecs := CanonicalGenerators(RootSystem(L))[3];
    BL := Basis(L);
    pos_slots := [];
    neg_slots := [];
    cartan_slots := [];
    for i in [1..Length(pos_vecs)] do
        Add(pos_slots, FindBasisSlotForLieVectorF4(BL, pos_vecs[i]));
    od;
    for i in [1..Length(neg_vecs)] do
        Add(neg_slots, FindBasisSlotForLieVectorF4(BL, neg_vecs[i]));
    od;
    for i in [1..Length(cartan_vecs)] do
        Add(cartan_slots, FindBasisSlotForLieVectorF4(BL, cartan_vecs[i]));
    od;
    if fail in pos_slots or fail in neg_slots or fail in cartan_slots then
        return fail;
    fi;
    return rec(
        basis := BL,
        positive_roots := pos_roots,
        positive_slots := pos_slots,
        negative_slots := neg_slots,
        cartan_slots := cartan_slots
    );
end;;

ReorderRank2CRootsForThetaF4 := function(subsystem_roots)
    local parser, root_vecs, lens;
    if not (IsList(subsystem_roots) and Length(subsystem_roots) = 2) then
        return fail;
    fi;
    if IsBoundGlobal("ParseRootStringFlexible") then
        parser := ValueGlobal("ParseRootStringFlexible");
    elif IsBoundGlobal("ParseRootString") then
        parser := ValueGlobal("ParseRootString");
    else
        return fail;
    fi;
    root_vecs := List(subsystem_roots, r -> parser(r));
    if fail in root_vecs then
        return fail;
    fi;
    lens := List(root_vecs, v -> F4InnerProductBySym(v, v));
    if lens[1] = fail or lens[2] = fail or lens[1] = lens[2] then
        return fail;
    fi;
    if lens[1] < lens[2] then
        return ShallowCopy(subsystem_roots);
    fi;
    return [subsystem_roots[2], subsystem_roots[1]];
end;;

MapSourceRootToSubsystemRootF4 := function(root_coeffs, subsystem_root_vecs)
    local out, i;
    if not IsList(root_coeffs) or not IsList(subsystem_root_vecs) or Length(root_coeffs) <> Length(subsystem_root_vecs) then
        return fail;
    fi;
    if Length(subsystem_root_vecs) = 0 then
        return fail;
    fi;
    out := List([1..Length(subsystem_root_vecs[1])], k -> 0);
    for i in [1..Length(root_coeffs)] do
        out := List([1..Length(out)], k -> out[k] + root_coeffs[i] * subsystem_root_vecs[i][k]);
    od;
    if IsBoundGlobal("IsF4RootVec") and not ValueGlobal("IsF4RootVec")(out) then
        return fail;
    fi;
    return ValueGlobal("FormatRootF4")(out);
end;;

AddThetaRootContributionF4 := function(e_map, f_map, h_map, source_roots, root, ex, fy)
    local hh;
    if ex = fail or fy = fail then
        return fail;
    fi;
    if ex <> 0 then
        Add(e_map, [root, ex]);
        Add(source_roots, root);
    fi;
    if fy <> 0 then
        Add(f_map, [root, fy]);
    fi;
    hh := ex * fy;
    if hh <> 0 then
        Add(h_map, [root, hh]);
    fi;
    return true;
end;;

AddThetaRootContributionFromSlotsF4 := function(x_coeffs, y_coeffs, x_slot, y_slot, root, e_map, f_map, h_map, source_roots)
    local ex, fy;
    ex := GetScaledCoefficientFromSlotF4(x_coeffs, x_slot);
    fy := GetScaledCoefficientFromSlotF4(y_coeffs, y_slot);
    return AddThetaRootContributionF4(e_map, f_map, h_map, source_roots, root, ex, fy);
end;;

BuildDirectThetaTripleRecordF4 := function(selected, e_map, f_map, h_map, source_roots)
    local mapped_h_expr, mapped_e_expr, mapped_f_expr;
    mapped_h_expr := FormatExprFromMapOrZeroF4(h_map, "H");
    mapped_e_expr := FormatExprFromMapOrZeroF4(e_map, "E");
    mapped_f_expr := FormatExprFromMapOrZeroF4(f_map, "F");
    if mapped_h_expr = "0" and mapped_e_expr = "0" and mapped_f_expr = "0" then
        return fail;
    fi;
    return rec(
        h := mapped_h_expr,
        e := mapped_e_expr,
        f := mapped_f_expr,
        source := "theta_noticed",
        source_x := String(selected.x),
        source_y := String(selected.y),
        source_theta_h := String(selected.h),
        source_roots := source_roots,
        source_direct_h := mapped_h_expr,
        source_direct_e := mapped_e_expr,
        source_direct_f := mapped_f_expr
    );
end;;

BuildDirectThetaTripleFromSourceDataF4 := function(selected, source_root_data, subsystem_root_vecs)
    local BL, x_coeffs, y_coeffs, e_map, f_map, h_map, source_roots, i, root_coeffs, mapped_root;
    BL := source_root_data.basis;
    x_coeffs := Coefficients(BL, selected.x);
    y_coeffs := Coefficients(BL, selected.y);
    e_map := [];
    f_map := [];
    h_map := [];
    source_roots := [];
    for i in [1..Length(source_root_data.positive_roots)] do
        root_coeffs := source_root_data.positive_roots[i];
        mapped_root := MapSourceRootToSubsystemRootF4(root_coeffs, subsystem_root_vecs);
        if mapped_root = fail then
            return fail;
        fi;
        if AddThetaRootContributionFromSlotsF4(
            x_coeffs, y_coeffs,
            source_root_data.positive_slots[i], source_root_data.negative_slots[i],
            mapped_root, e_map, f_map, h_map, source_roots
        ) = fail then
            return fail;
        fi;
        mapped_root := MapSourceRootToSubsystemRootF4(List(root_coeffs, c -> -c), subsystem_root_vecs);
        if mapped_root = fail then
            return fail;
        fi;
        if AddThetaRootContributionFromSlotsF4(
            x_coeffs, y_coeffs,
            source_root_data.negative_slots[i], source_root_data.positive_slots[i],
            mapped_root, e_map, f_map, h_map, source_roots
        ) = fail then
            return fail;
        fi;
    od;
    return BuildDirectThetaTripleRecordF4(selected, e_map, f_map, h_map, source_roots);
end;;

BuildDirectThetaTripleForF4 := function(theta_data, subsystem_roots, colors_opt)
    local selected, real_form, BL, x_coeffs, y_coeffs, direct_data, slot_map, slot, e_map, f_map, h_map,
          parser, subsystem_root_vecs, source_root_data, root_coeffs, mapped_root, i, source_roots,
          reordered_roots, reordered_root_vecs;
    if theta_data = fail or not IsRecord(theta_data) or not IsBound(theta_data.ok) or not theta_data.ok then
        return fail;
    fi;
    if not IsBound(theta_data.selected) or theta_data.selected = fail then
        return fail;
    fi;
    if not IsBound(theta_data.real_form) or theta_data.real_form = fail then
        return fail;
    fi;
    if not IsBound(theta_data.theta) or theta_data.theta = fail then
        return fail;
    fi;
    selected := theta_data.selected;
    real_form := theta_data.real_form;
    source_root_data := fail;
    if real_form.type = "sl_R" and Length(subsystem_roots) >= 1 then
        source_root_data := BuildThetaSourceRootDataF4(theta_data);
        if source_root_data <> fail and Length(source_root_data.cartan_slots) = 1 then
            BL := source_root_data.basis;
            x_coeffs := Coefficients(BL, selected.x);
            y_coeffs := Coefficients(BL, selected.y);
            e_map := [];
            f_map := [];
            h_map := [];
            source_roots := [];
            if AddThetaRootContributionFromSlotsF4(
                x_coeffs, y_coeffs,
                source_root_data.positive_slots[1], source_root_data.negative_slots[1],
                subsystem_roots[1], e_map, f_map, h_map, source_roots
            ) = fail then
                return fail;
            fi;
            if AddThetaRootContributionFromSlotsF4(
                x_coeffs, y_coeffs,
                source_root_data.negative_slots[1], source_root_data.positive_slots[1],
                NegateRootStringF4(subsystem_roots[1]), e_map, f_map, h_map, source_roots
            ) = fail then
                return fail;
            fi;
            return BuildDirectThetaTripleRecordF4(selected, e_map, f_map, h_map, source_roots);
        fi;
        BL := Basis(Source(theta_data.theta));
        x_coeffs := Coefficients(BL, selected.x);
        y_coeffs := Coefficients(BL, selected.y);
        slot_map := [
            rec(root := subsystem_roots[1], x_index := 2, y_index := 1),
            rec(root := NegateRootStringF4(subsystem_roots[1]), x_index := 1, y_index := 2)
        ];
        e_map := [];
        f_map := [];
        h_map := [];
        source_roots := [];
        for slot in slot_map do
            if AddThetaRootContributionF4(
                e_map, f_map, h_map, source_roots,
                slot.root, x_coeffs[slot.x_index], y_coeffs[slot.y_index]
            ) = fail then
                return fail;
            fi;
        od;
        return BuildDirectThetaTripleRecordF4(selected, e_map, f_map, h_map, source_roots);
    elif (
        real_form.type = "su_pq" or real_form.type = "so_pq" or
        real_form.type = "sp_R" or real_form.type = "sp_pq"
    ) and IsBound(real_form.params) and Length(real_form.params) >= 1 then
        if not (IsList(subsystem_roots) and Length(subsystem_roots) > 0) then
            return fail;
        fi;
        if IsBoundGlobal("ParseRootStringFlexible") then
            parser := ValueGlobal("ParseRootStringFlexible");
        elif IsBoundGlobal("ParseRootString") then
            parser := ValueGlobal("ParseRootString");
        else
            return fail;
        fi;
        subsystem_root_vecs := List(subsystem_roots, r -> parser(r));
        if fail in subsystem_root_vecs then
            return fail;
        fi;
        source_root_data := BuildThetaSourceRootDataF4(theta_data);
        if source_root_data = fail then
            if not (real_form.type = "su_pq" and real_form.params[1] = 2 and real_form.params[2] = 2) then
                return fail;
            fi;
            direct_data := BuildDirectThetaRootListForF4(theta_data, subsystem_roots, colors_opt);
            if direct_data = fail or not IsRecord(direct_data) or not IsBound(direct_data.roots) or Length(direct_data.roots) <> 3 then
                return fail;
            fi;
            BL := Basis(Source(theta_data.theta));
            x_coeffs := Coefficients(BL, selected.x);
            y_coeffs := Coefficients(BL, selected.y);
            if IsBound(direct_data.negated) and direct_data.negated then
                slot_map := [
                    rec(root := direct_data.roots[1], x_index := 6, y_index := 4),
                    rec(root := direct_data.roots[2], x_index := 10, y_index := 5),
                    rec(root := direct_data.roots[3], x_index := 11, y_index := 12)
                ];
            else
                slot_map := [
                    rec(root := direct_data.roots[1], x_index := 4, y_index := 2),
                    rec(root := direct_data.roots[2], x_index := 5, y_index := 10),
                    rec(root := direct_data.roots[3], x_index := 8, y_index := 11)
                ];
            fi;
            e_map := [];
            f_map := [];
            h_map := [];
            source_roots := [];
            for slot in slot_map do
                if AddThetaRootContributionF4(
                    e_map, f_map, h_map, source_roots,
                    slot.root, x_coeffs[slot.x_index], y_coeffs[slot.y_index]
                ) = fail then
                    return fail;
                fi;
            od;
            return BuildDirectThetaTripleRecordF4(selected, e_map, f_map, h_map, source_roots);
        fi;
        BL := source_root_data.basis;
        direct_data := BuildDirectThetaTripleFromSourceDataF4(selected, source_root_data, subsystem_root_vecs);
        if direct_data <> fail then
            return direct_data;
        fi;
        if Length(subsystem_roots) = 2 and (real_form.type = "sp_R" or real_form.type = "sp_pq") then
            reordered_roots := ReorderRank2CRootsForThetaF4(subsystem_roots);
            if reordered_roots <> fail and reordered_roots <> subsystem_roots then
                reordered_root_vecs := List(reordered_roots, r -> parser(r));
                if not (fail in reordered_root_vecs) then
                    direct_data := BuildDirectThetaTripleFromSourceDataF4(selected, source_root_data, reordered_root_vecs);
                    if direct_data <> fail then
                        return direct_data;
                    fi;
                fi;
            fi;
        fi;
        return fail;
    else
        return fail;
    fi;
end;;

TranslateSelectedThetaTripleToRootTripleF4 := function(theta_data, subsystem_roots, colors_opt)
    local direct_triple;
    if theta_data = fail or not IsRecord(theta_data) or not IsBound(theta_data.ok) or not theta_data.ok then
        return fail;
    fi;
    if not IsBound(theta_data.selected) or theta_data.selected = fail then
        return fail;
    fi;
    if not IsBound(theta_data.real_form) or theta_data.real_form = fail then
        return fail;
    fi;
    direct_triple := BuildDirectThetaTripleForF4(theta_data, subsystem_roots, colors_opt);
    if direct_triple <> fail then
        return direct_triple;
    fi;
    return fail;
end;;

CollectThetaRootTripleCandidatesF4 := function(theta_data, subsystem_roots, colors_opt)
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
                one := TranslateSelectedThetaTripleToRootTripleF4(td, subsystem_roots, colors_opt);
                if one <> fail then
                    Add(candidates, rec(triple := one, orbit_index := t.index));
                fi;
            fi;
        od;
    fi;
    return candidates;
end;;

ReorderCompositeComponentsForDisplay := function(component_data)
    local indexed, i, item, j;
    indexed := [];
    for i in [1..Length(component_data)] do
        Add(indexed, rec(idx := i, data := component_data[i], rank := Length(component_data[i].roots)));
    od;
    for i in [2..Length(indexed)] do
        item := indexed[i];
        j := i - 1;
        while j >= 1 and indexed[j].rank < item.rank do
            indexed[j + 1] := indexed[j];
            j := j - 1;
        od;
        indexed[j + 1] := item;
    od;
    return List(indexed, x -> x.data);
end;;

VerifyThetaSelectedTripleOnLieLevelF4 := function(theta_data)
    local selected, hx_ok, hy_ok, xy_ok;
    if theta_data = fail or not IsRecord(theta_data) then
        return rec(ok := false, reason := "theta_data_fail");
    fi;
    if not IsBound(theta_data.selected) or theta_data.selected = fail or not IsRecord(theta_data.selected) then
        return rec(ok := false, reason := "selected_fail");
    fi;
    selected := theta_data.selected;
    if not IsBound(selected.x) or not IsBound(selected.y) or not IsBound(selected.h) then
        return rec(ok := false, reason := "selected_xyz_missing");
    fi;
    hx_ok := IsZero(selected.h * selected.x - 2 * selected.x);
    hy_ok := IsZero(selected.h * selected.y + 2 * selected.y);
    xy_ok := IsZero(selected.x * selected.y - selected.h);
    return rec(ok := (hx_ok and hy_ok and xy_ok), hx := hx_ok, hy := hy_ok, xy := xy_ok);
end;;

PrintSL2TripleUnified := function(type, rank, subsystem_roots, colors_opt)
    local normalized_roots, e_str, h_str, f_str, is_bdi_type, is_aiii_type, effective_colors, triple, ambient, efh_ok, pk_ok, c1, c2, c3, c4, f4_fix, diag_ok, real_form_info, theta_data, theta_triple, is_composite, k, m5_item, m5_count, theta_obj_check, theta_candidates, theta_candidate_idx, theta_candidate_item, raw_m5_count;

    normalized_roots := NormalizeRootsForNormalTriple(type, rank, subsystem_roots, colors_opt);

    

    is_bdi_type := (type = "BDI" or type = "so_pq");
    is_aiii_type := (type = "su_pq" or type = "AIII" or PositionSublist(type, "su(") <> fail or PositionSublist(type, "AIII") <> fail);
    effective_colors := colors_opt;
    if (not is_bdi_type) then
        if rank = 4 and IsList(effective_colors) and Length(effective_colors) >= 4 then
            is_bdi_type := true;
        fi;
    fi;
    if rank = 4 and IsList(effective_colors) and Length(effective_colors) = 4 and Position(subsystem_roots, "-2a1-3a2-4a3-2a4") <> fail then
        if IsBoundGlobal("IsNoncompactF4Vec") then
            c1 := 0; c2 := 0; c3 := 0; c4 := 0;
            if ValueGlobal("IsNoncompactF4Vec")([1,0,0,0]) then c1 := 1; fi;
            if ValueGlobal("IsNoncompactF4Vec")([0,1,0,0]) then c2 := 1; fi;
            if ValueGlobal("IsNoncompactF4Vec")([0,0,1,0]) then c3 := 1; fi;
            if ValueGlobal("IsNoncompactF4Vec")([0,0,0,1]) then c4 := 1; fi;
            effective_colors := [c1, c2, c3, c4];
        fi;
    fi;
    if (not IsList(effective_colors) or Length(effective_colors) = 0) and IsBoundGlobal("GlobalColorsListForDiagram") then
        effective_colors := ValueGlobal("GlobalColorsListForDiagram");
    fi;
    if (not is_bdi_type) then
        if rank = 4 and IsList(effective_colors) and Length(effective_colors) >= 4 then
            is_bdi_type := true;
        fi;
    fi;

    ambient := "F4";
    if IsBoundGlobal("GlobalAmbientType") then
        ambient := ValueGlobal("GlobalAmbientType");
    fi;
    if rank = 4 and (type = "B" or type = "BDI") then
        if IsBoundGlobal("SanitizeF4RootList") then
            normalized_roots := ValueGlobal("SanitizeF4RootList")(normalized_roots);
        fi;
    fi;

    is_composite := (type = "Composite" and IsBoundGlobal("GlobalCompositeComponentData") and IsList(ValueGlobal("GlobalCompositeComponentData")) and Length(ValueGlobal("GlobalCompositeComponentData")) > 0);
    theta_data := fail;
    theta_triple := fail;
    theta_candidates := [];
    if is_composite and IsBoundGlobal("BuildCompositeNormalTriple_FromComponents") then
        GlobalCurrentThetaOrbitData := fail;
        triple := ValueGlobal("BuildCompositeNormalTriple_FromComponents")(ValueGlobal("GlobalCompositeComponentData"));
        if IsBoundGlobal("GlobalF4CompositeOutputStyle") and ValueGlobal("GlobalF4CompositeOutputStyle") = "single_like" then
            Print("        [M4-STEP3] noticed triple 已按组件合成到根形式\n");
        else
            Print("        [组件求和-RAW] h = ", triple.h, "\n");
            Print("        [组件求和-RAW] e = ", triple.e, "\n");
            Print("        [组件求和-RAW] f = ", triple.f, "\n");
            Print("        [组件求和] 按组件构造并直接求和完成\n");
        fi;
    elif ambient = "F4" then
        real_form_info := ResolveModule4RealFormInfo(type, rank);
        theta_data := ComputeThetaOrbitTriplesForModule4(real_form_info);
        GlobalCurrentThetaOrbitData := theta_data;
        PrintThetaOrbitStep2(theta_data);
        if theta_data <> fail and IsRecord(theta_data) and IsBound(theta_data.ok) and theta_data.ok then
            theta_candidates := CollectThetaRootTripleCandidatesF4(theta_data, normalized_roots, effective_colors);
            theta_triple := TranslateSelectedThetaTripleToRootTripleF4(theta_data, normalized_roots, effective_colors);
        fi;
        if Length(theta_candidates) > 0 then
            Print("        [M4-STEP3] noticed 候选数 = ", Length(theta_candidates), "\n");
            for theta_candidate_idx in [1..Length(theta_candidates)] do
                theta_candidate_item := theta_candidates[theta_candidate_idx];
                Print("        [THETA-CANDIDATE#", theta_candidate_idx, "] h = ", theta_candidate_item.triple.h, "\n");
                Print("        [THETA-CANDIDATE#", theta_candidate_idx, "] e = ", theta_candidate_item.triple.e, "\n");
                Print("        [THETA-CANDIDATE#", theta_candidate_idx, "] f = ", theta_candidate_item.triple.f, "\n");
                if IsBound(theta_candidate_item.orbit_index) then
                    Print("        [THETA-CANDIDATE#", theta_candidate_idx, "] 来自 orbit#", theta_candidate_item.orbit_index, "\n");
                fi;
            od;
        fi;
        if theta_triple <> fail then
            triple := theta_triple;
            if Length(theta_candidates) > 0 then
                triple.all_triples := List(theta_candidates, x -> rec(h := x.triple.h, e := x.triple.e, f := x.triple.f, orbit_index := x.orbit_index));
            fi;
            Print("        [M4-STEP3] noticed triple 已翻译到根形式\n");
            if IsBound(triple.source_x) and IsBound(triple.source_roots) then
                Print("        [M4-STEP3] 根位匹配: x=", triple.source_x, " => roots=", triple.source_roots, "\n");
            fi;
            if IsBound(triple.source_x) and IsBound(triple.source_direct_e) then
                Print("        [M4-STEP3] 候选根式: x=", triple.source_x, " => e=", triple.source_direct_e, "\n");
            fi;
            if IsBound(triple.source_y) and IsBound(triple.source_direct_f) then
                Print("        [M4-STEP3] 候选根式: y=", triple.source_y, " => f=", triple.source_direct_f, "\n");
            fi;
            if IsBound(triple.source_theta_h) and IsBound(triple.source_direct_h) then
                Print("        [M4-STEP3] 由e,f重建: h=", triple.source_theta_h, " => H=", triple.source_direct_h, "\n");
            fi;
            Print("        [THETA-CANDIDATE] h = ", triple.h, "\n");
            Print("        [THETA-CANDIDATE] e = ", triple.e, "\n");
            Print("        [THETA-CANDIDATE] f = ", triple.f, "\n");
        else
            Error("F4 主路径要求 theta noticed triple 必须成功翻译到根形式。");
        fi;
    else
        Error("F4 主流程仅保留 theta-orbit 构造路径。");
    fi;
    h_str := triple.h;
    e_str := triple.e;
    f_str := triple.f;
    Print("        >>> 模块 4 输出: SL2-Triple 构造 (Candidate)\n");
    Print("        SL2-Triple (h, e, f):\n");
    Print("          h = ", h_str, "\n");
    Print("          e = ", e_str, "\n");
    Print("          f = ", f_str, "\n");

    Print("        [强制校验] 开始 VerifyEFH_NormalConditions…\n");
    if IsBoundGlobal("VerifyEFH_NormalConditions") then
        efh_ok := ValueGlobal("VerifyEFH_NormalConditions")(h_str, e_str, f_str);
    else
        efh_ok := false;
        Print("        [强制校验] VerifyEFH_NormalConditions 未加载！\n");
    fi;
    if ambient = "F4" and IsRecord(triple) and IsBound(triple.source) and triple.source = "theta_noticed" then
        theta_obj_check := VerifyThetaSelectedTripleOnLieLevelF4(theta_data);
        if IsRecord(theta_obj_check) and IsBound(theta_obj_check.hx) and IsBound(theta_obj_check.hy) and IsBound(theta_obj_check.xy) then
            Print("        [Theta对象校验] [h,x]=2x: ", theta_obj_check.hx, "\n");
            Print("        [Theta对象校验] [h,y]=-2y: ", theta_obj_check.hy, "\n");
            Print("        [Theta对象校验] [x,y]=h: ", theta_obj_check.xy, "\n");
        fi;
        if IsRecord(theta_obj_check) and IsBound(theta_obj_check.ok) then
            efh_ok := theta_obj_check.ok;
            Print("        [强制校验] 模块4采用 Theta 对象校验结果\n");
        fi;
    fi;
    pk_ok := true;
    if ambient = "G2" and IsBoundGlobal("VerifyPK_Subsystem_G2") then
        pk_ok := ValueGlobal("VerifyPK_Subsystem_G2")(h_str, e_str, f_str, normalized_roots, effective_colors);
        if pk_ok then
            Print("        [诊断] P/K 检查: PASS\n");
        else
            Print("        [诊断] P/K 检查: FAIL\n");
        fi;
    elif ambient = "F4" and IsBoundGlobal("VerifyPK_Subsystem_F4") then
        pk_ok := ValueGlobal("VerifyPK_Subsystem_F4")(h_str, e_str, f_str, normalized_roots, effective_colors);
        if pk_ok then
            Print("        [诊断] P/K 检查: PASS\n");
        else
            Print("        [诊断] P/K 检查: FAIL\n");
        fi;
    fi;
    if ambient = "F4" and type = "Composite" and (not efh_ok or not pk_ok) then
        f4_fix := RefitCompositeTripleF4_PostCheck(h_str, e_str, f_str, normalized_roots, effective_colors, efh_ok, pk_ok);
        h_str := f4_fix.h;
        e_str := f4_fix.e;
        f_str := f4_fix.f;
        efh_ok := f4_fix.efh_ok;
        pk_ok := f4_fix.pk_ok;
    fi;
    if ambient = "F4" then
        diag_ok := efh_ok;
    else
        diag_ok := (efh_ok and pk_ok);
    fi;
    if diag_ok then
        Print("        模块4判定: Normal\n");
        Print("NORMAL_TRIPLE_CHECK: PASS\n");
        RecordNormalTripleStats(true);
    else
        Print("        模块4判定: Non-normal\n");
        Print("NORMAL_TRIPLE_CHECK: FAIL\n");
        RecordNormalTripleStats(false);
    fi;

    Print("\n");
    if not IsBoundGlobal("PrintKOrbitInfo_F4") and not IsBoundGlobal("PrintKOrbitInfo") and IsBoundGlobal("Read") then
        Read("../03_Orbit_Analysis/K_Orbit_Classifier.g");
    fi;
    ambient := "F4";
    if IsBoundGlobal("GlobalAmbientType") then
        ambient := ValueGlobal("GlobalAmbientType");
    fi;
    if IsRecord(triple) and IsBound(triple.all_triples) and IsList(triple.all_triples) and Length(triple.all_triples) > 0 then
        raw_m5_count := Length(triple.all_triples);
        if ambient = "F4" and IsBoundGlobal("FilterTriplesByWKConjugacyFixingDeltaJF4") then
            triple.all_triples := ValueGlobal("FilterTriplesByWKConjugacyFixingDeltaJF4")(triple.all_triples, subsystem_roots);
        fi;
        if ambient = "F4" and IsBoundGlobal("FilterCompositeTripleFamiliesByKOrbitLabelF4") then
            triple.all_triples := ValueGlobal("FilterCompositeTripleFamiliesByKOrbitLabelF4")(triple.all_triples);
        fi;
        m5_count := Length(triple.all_triples);
        if ambient = "F4" then
            if m5_count < raw_m5_count then
                Print("        [M4-WK-FILTER] 过滤后 noticed triple 数: ", raw_m5_count, " -> ", m5_count, "\n");
            else
                Print("        [M4-WK-FILTER] 无可过滤的 noticed triple\n");
            fi;
        fi;
        if is_composite then
            if IsBoundGlobal("GlobalF4CompositeOutputStyle") and ValueGlobal("GlobalF4CompositeOutputStyle") = "single_like" then
                Print("        [M5-PROPAGATE] 向模块5传递 noticed triple 数 = ", m5_count, "\n");
            else
                Print("        [M5-PROPAGATE] 向模块5传递组合 triple 数 = ", m5_count, "\n");
            fi;
        else
            Print("        [M5-PROPAGATE] 向模块5传递 noticed triple 数 = ", m5_count, "\n");
        fi;
        for k in [1..m5_count] do
            m5_item := triple.all_triples[k];
            Print("        [M5-PROPAGATE#", k, "] h = ", m5_item.h, "\n");
            Print("        [M5-PROPAGATE#", k, "] e = ", m5_item.e, "\n");
            Print("        [M5-PROPAGATE#", k, "] f = ", m5_item.f, "\n");
            if IsBound(m5_item.component_orbit_choice) then
                Print("        [M5-PROPAGATE#", k, "] 组件orbit选择 = ", m5_item.component_orbit_choice, "\n");
            fi;
            if IsBound(m5_item.orbit_index) then
                Print("        [M5-PROPAGATE#", k, "] 来源 orbit#", m5_item.orbit_index, "\n");
            fi;
            DispatchKOrbitOutput(ambient, type, m5_item.h, m5_item.e, m5_item.f);
        od;
    else
        DispatchKOrbitOutput(ambient, type, h_str, e_str, f_str);
    fi;
end;;

# =============================================================================
# 核心构造逻辑
# =============================================================================

# ConstructEElement: 构造幂零元 e
# 当前主流程已改为 theta-orbit 直译，这里保留的是根名称映射说明。

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
# 当前主流程由 theta-orbit 给出 h,e,f，这里保留的是旧链式解释的背景注释。
# =============================================================================
# 主接口：打印 (h, e, f) 并联动 K-Orbit 分类输出
# =============================================================================

PrintGenericHSummary := function(h_expr)
    local parts, i, t, coeffs, c, ch, num, sign, k, term, coef, inside, ambient;
    ambient := "F4";
    if IsBoundGlobal("GlobalAmbientType") then
        ambient := ValueGlobal("GlobalAmbientType");
    fi;
    if ambient = "G2" and IsBoundGlobal("PrintKOrbitInfo") then
        GlobalCurrentSubsystemType := "Unknown";
        ValueGlobal("PrintKOrbitInfo")(h_expr);
        return;
    elif IsBoundGlobal("PrintKOrbitInfo_F4") then
        ValueGlobal("PrintKOrbitInfo_F4")(h_expr, [0,0,0,0,0]);
        return;
    elif IsBoundGlobal("DrawF4ExtendedVoganDiagram") then
        Print("        >>> 模块 5 输出: K-Orbit 分类\n");
        Print("        K-Orbit Classification (F4(4)):\n");
        Print("          初始 h 系数 (F4 扩展): <auto-parse>\n");
        ValueGlobal("DrawF4ExtendedVoganDiagram")(h_expr, [0,0,0,0,0]);
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

BuildCompositeNormalTriple_FromComponents := function(component_data)
    local comp, comp_roots, comp_colors, comp_type, comp_rank, comp_triple,
          h_expr, e_expr, f_expr, normalized_comp_roots, has_su_pq,
          ambient, comp_efh_ok, comp_pk_ok, comp_diag_ok, comp_fix,
          comp_real_form_info, theta_data, comp_candidates, candidate, cand_idx,
          all_component_candidates, combined, new_combined, base, item, merged, path, k, comp_idx,
          ordered_component_data, composite_single_like;
    h_expr := "";
    e_expr := "";
    f_expr := "";
    has_su_pq := false;
    all_component_candidates := [];
    composite_single_like := (IsBoundGlobal("GlobalF4CompositeOutputStyle") and ValueGlobal("GlobalF4CompositeOutputStyle") = "single_like");
    ordered_component_data := ReorderCompositeComponentsForDisplay(component_data);
    comp_idx := 0;
    for comp in ordered_component_data do
        comp_idx := comp_idx + 1;
        comp_roots := comp.roots;
        comp_colors := fail;
        if IsBound(comp.colors) then
            comp_colors := comp.colors;
        fi;
        comp_type := "CompositeComponent";
        if IsBound(comp.real_form_type) then
            comp_type := comp.real_form_type;
        fi;
        if comp_type = "su_pq" or PositionSublist(comp_type, "su(") <> fail or PositionSublist(comp_type, "AIII") <> fail then
            has_su_pq := true;
        fi;
        comp_rank := Length(comp_roots);
        normalized_comp_roots := NormalizeRootsForNormalTriple(comp_type, comp_rank, comp_roots, comp_colors);
        if IsBoundGlobal("GlobalAmbientType") and ValueGlobal("GlobalAmbientType") = "G2" and (comp_type = "su_pq" or PositionSublist(comp_type, "su(") <> fail or PositionSublist(comp_type, "AIII") <> fail) then
            normalized_comp_roots := comp_roots;
        fi;
        ambient := "F4";
        if IsBoundGlobal("GlobalAmbientType") then
            ambient := ValueGlobal("GlobalAmbientType");
        fi;
        if not composite_single_like then
            Print("        [COMP#", comp_idx, "] 开始模块4组件流程\n");
            Print("        [组件诊断] 使用theta-orbit构造\n");
        fi;
        comp_real_form_info := fail;
        if IsBound(comp.real_form_type) and IsBound(comp.params) then
            comp_real_form_info := rec(type := comp.real_form_type, params := comp.params);
        fi;
        theta_data := ComputeThetaOrbitTriplesForModule4(comp_real_form_info);
        GlobalCurrentThetaOrbitData := theta_data;
        if not composite_single_like then
            PrintThetaOrbitStep2(theta_data);
        fi;
        comp_candidates := CollectThetaRootTripleCandidatesF4(theta_data, normalized_comp_roots, comp_colors);
        if Length(comp_candidates) = 0 then
            Error("F4 复合主路径要求每个组件的 theta noticed triple 都能成功翻译到根形式。");
        fi;
        if not composite_single_like then
            Print("        [组件求和输入] 来源: noticed triple 根形式翻译\n");
            Print("        [COMP#", comp_idx, "] noticed 候选数 = ", Length(comp_candidates), "\n");
        fi;
        cand_idx := 0;
        for candidate in comp_candidates do
            cand_idx := cand_idx + 1;
            if not composite_single_like then
                Print("        [组件候选#", cand_idx, "] h = ", candidate.triple.h, "\n");
                Print("        [组件候选#", cand_idx, "] e = ", candidate.triple.e, "\n");
                Print("        [组件候选#", cand_idx, "] f = ", candidate.triple.f, "\n");
                if IsBound(candidate.triple.source_x) and IsBound(candidate.triple.source_roots) then
                    Print("        [组件候选#", cand_idx, "] 根位匹配: x=", candidate.triple.source_x, " => roots=", candidate.triple.source_roots, "\n");
                fi;
                if IsBound(candidate.triple.source_x) and IsBound(candidate.triple.source_direct_e) then
                    Print("        [组件候选#", cand_idx, "] 候选根式: x=", candidate.triple.source_x, " => e=", candidate.triple.source_direct_e, "\n");
                fi;
                if IsBound(candidate.triple.source_y) and IsBound(candidate.triple.source_direct_f) then
                    Print("        [组件候选#", cand_idx, "] 候选根式: y=", candidate.triple.source_y, " => f=", candidate.triple.source_direct_f, "\n");
                fi;
                if IsBound(candidate.triple.source_theta_h) and IsBound(candidate.triple.source_direct_h) then
                    Print("        [组件候选#", cand_idx, "] 由e,f重建: h=", candidate.triple.source_theta_h, " => H=", candidate.triple.source_direct_h, "\n");
                fi;
                if candidate.orbit_index <> fail then
                    Print("        [组件候选#", cand_idx, "] 来自 orbit#", candidate.orbit_index, "\n");
                fi;
            fi;
        od;
        comp_triple := comp_candidates[1].triple;
        if not composite_single_like then
            Print("        [组件诊断] 类型=", comp_type, ", 根=", comp_roots, "\n");
            Print("        [组件诊断] 开始 VerifyEFH_NormalConditions…\n");
        fi;
        comp_efh_ok := false;
        if IsBoundGlobal("VerifyEFH_NormalConditions") then
            comp_efh_ok := ValueGlobal("VerifyEFH_NormalConditions")(comp_triple.h, comp_triple.e, comp_triple.f);
        fi;
        comp_pk_ok := true;
        if ambient = "G2" and IsBoundGlobal("VerifyPK_Subsystem_G2") then
            comp_pk_ok := ValueGlobal("VerifyPK_Subsystem_G2")(comp_triple.h, comp_triple.e, comp_triple.f, normalized_comp_roots, comp_colors);
            if not composite_single_like then
                if comp_pk_ok then
                    Print("        [组件诊断] P/K 检查: PASS\n");
                else
                    Print("        [组件诊断] P/K 检查: FAIL\n");
                fi;
            fi;
        elif ambient = "F4" and IsBoundGlobal("VerifyPK_Subsystem_F4") then
            comp_pk_ok := ValueGlobal("VerifyPK_Subsystem_F4")(comp_triple.h, comp_triple.e, comp_triple.f, normalized_comp_roots, comp_colors);
            if not composite_single_like then
                if comp_pk_ok then
                    Print("        [组件诊断] P/K 检查: PASS\n");
                else
                    Print("        [组件诊断] P/K 检查: FAIL\n");
                fi;
            fi;
        fi;
        if ambient = "F4" and (not comp_efh_ok or not comp_pk_ok) then
            comp_fix := RefitCompositeTripleF4_PostCheck(comp_triple.h, comp_triple.e, comp_triple.f, normalized_comp_roots, comp_colors, comp_efh_ok, comp_pk_ok);
            comp_triple := rec(h := comp_fix.h, e := comp_fix.e, f := comp_fix.f);
            comp_efh_ok := comp_fix.efh_ok;
            comp_pk_ok := comp_fix.pk_ok;
        fi;
        if ambient = "F4" then
            comp_diag_ok := comp_efh_ok;
        else
            comp_diag_ok := (comp_efh_ok and comp_pk_ok);
        fi;
        if not composite_single_like then
            if comp_diag_ok then
                Print("        [组件诊断] NORMAL_TRIPLE_CHECK: PASS\n");
            else
                Print("        [组件诊断] NORMAL_TRIPLE_CHECK: FAIL\n");
            fi;
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
    if composite_single_like then
        Print("        [M4-STEP2] noticed 轨道数 = ", Length(combined), "\n");
    else
        Print("        [COMPOSITE] 组件候选笛卡尔组合数 = ", Length(combined), "\n");
    fi;
    for k in [1..Length(combined)] do
        if combined[k].h = "" then combined[k].h := "0"; fi;
        if combined[k].e = "" then combined[k].e := "0"; fi;
        if combined[k].f = "" then combined[k].f := "0"; fi;
        if composite_single_like then
            Print("        [THETA-CANDIDATE#", k, "] h = ", combined[k].h, "\n");
            Print("        [THETA-CANDIDATE#", k, "] e = ", combined[k].e, "\n");
            Print("        [THETA-CANDIDATE#", k, "] f = ", combined[k].f, "\n");
            Print("        [THETA-CANDIDATE#", k, "] 组件orbit选择 = ", combined[k].path, "\n");
        else
            Print("        [COMPOSITE-RAW#", k, "] h = ", combined[k].h, "\n");
            Print("        [COMPOSITE-RAW#", k, "] e = ", combined[k].e, "\n");
            Print("        [COMPOSITE-RAW#", k, "] f = ", combined[k].f, "\n");
            Print("        [COMPOSITE-RAW#", k, "] 组件orbit选择 = ", combined[k].path, "\n");
        fi;
    od;
    h_expr := combined[1].h;
    e_expr := combined[1].e;
    f_expr := combined[1].f;
    if h_expr = "" then h_expr := "0"; fi;
    if e_expr = "" then e_expr := "0"; fi;
    if f_expr = "" then f_expr := "0"; fi;
    if has_su_pq and IsBoundGlobal("GlobalAmbientType") and ValueGlobal("GlobalAmbientType") = "F4" then
        Print("        [组件求和] AIII 全局重整: SKIPPED (RAW=各组件直接相加)\n");
    fi;
    return rec(h := h_expr, e := e_expr, f := f_expr, all_triples := List(combined, c -> rec(h := c.h, e := c.e, f := c.f, component_orbit_choice := c.path)));
end;;

ZeroSquareMatrix := function(n)
    return List([1..n], i -> List([1..n], j -> 0));
end;;

DiagonalMatrixFromList := function(diag)
    local n, i, m;
    n := Length(diag);
    m := ZeroSquareMatrix(n);
    for i in [1..n] do
        m[i][i] := diag[i];
    od;
    return m;
end;;

MatrixCommutator := function(x, y)
    return x * y - y * x;
end;;

IsInFAIII := function(x, dim_va, dim_vb)
    local i, j, n;
    n := dim_va + dim_vb;
    if Length(x) <> n then
        return false;
    fi;
    for i in [1..n] do
        if Length(x[i]) <> n then
            return false;
        fi;
    od;
    for i in [1..dim_va] do
        for j in [dim_va + 1..n] do
            if x[i][j] <> 0 then
                return false;
            fi;
        od;
    od;
    for i in [dim_va + 1..n] do
        for j in [1..dim_va] do
            if x[i][j] <> 0 then
                return false;
            fi;
        od;
    od;
    return true;
end;;

IsInPAIII := function(x, dim_va, dim_vb)
    local i, j, n;
    n := dim_va + dim_vb;
    if Length(x) <> n then
        return false;
    fi;
    for i in [1..n] do
        if Length(x[i]) <> n then
            return false;
        fi;
    od;
    for i in [1..dim_va] do
        for j in [1..dim_va] do
            if x[i][j] <> 0 then
                return false;
            fi;
        od;
    od;
    for i in [dim_va + 1..n] do
        for j in [dim_va + 1..n] do
            if x[i][j] <> 0 then
                return false;
            fi;
        od;
    od;
    return true;
end;;

BuildAIIIWeightData := function(dim_va, dim_vb)
    local alpha_names, beta_names, lambda_va, lambda_vb, delta_g, delta_f, delta_p, i, j;
    alpha_names := [];
    beta_names := [];
    for i in [1..dim_va] do
        Add(alpha_names, Concatenation("alpha_", String(i)));
    od;
    for j in [1..dim_vb] do
        Add(beta_names, Concatenation("beta_", String(j)));
    od;
    lambda_va := ShallowCopy(alpha_names);
    lambda_vb := ShallowCopy(beta_names);
    delta_g := [];
    for i in [1..dim_va] do
        for j in [1..dim_va] do
            if i <> j then
                Add(delta_g, Concatenation(alpha_names[i], " - ", alpha_names[j]));
            fi;
        od;
    od;
    for i in [1..dim_vb] do
        for j in [1..dim_vb] do
            if i <> j then
                Add(delta_g, Concatenation(beta_names[i], " - ", beta_names[j]));
            fi;
        od;
    od;
    for i in [1..dim_va] do
        for j in [1..dim_vb] do
            Add(delta_g, Concatenation(alpha_names[i], " - ", beta_names[j]));
            Add(delta_g, Concatenation(beta_names[j], " - ", alpha_names[i]));
        od;
    od;
    delta_f := [];
    for i in [1..dim_va] do
        for j in [1..dim_va] do
            if i <> j then
                Add(delta_f, Concatenation(alpha_names[i], " - ", alpha_names[j]));
            fi;
        od;
    od;
    for i in [1..dim_vb] do
        for j in [1..dim_vb] do
            if i <> j then
                Add(delta_f, Concatenation(beta_names[i], " - ", beta_names[j]));
            fi;
        od;
    od;
    delta_p := [];
    for i in [1..dim_va] do
        for j in [1..dim_vb] do
            Add(delta_p, Concatenation(alpha_names[i], " - ", beta_names[j]));
            Add(delta_p, Concatenation(beta_names[j], " - ", alpha_names[i]));
        od;
    od;
    return rec(
        t_desc := "{t in g; acts on each a_i,b_j as a scalar}",
        alpha_names := alpha_names,
        beta_names := beta_names,
        alpha_action := List([1..dim_va], i -> Concatenation("t a_", String(i), " = alpha_", String(i), "(t) a_", String(i))),
        beta_action := List([1..dim_vb], j -> Concatenation("t b_", String(j), " = beta_", String(j), "(t) b_", String(j))),
        Lambda_Va_t := lambda_va,
        Lambda_Vb_t := lambda_vb,
        Delta_g_t := delta_g,
        Delta_f_t := delta_f,
        Delta_p_t := delta_p
    );
end;;

EvaluateAIIIWeightsOnH := function(h_mat, basis_sequence, dim_va, dim_vb)
    local n, i, diag_vals, alpha_vals, beta_vals, name, idx_str, idx, good, delta_f_vals, delta_p_vals, j;
    n := dim_va + dim_vb;
    if Length(h_mat) <> n then
        return rec(ok := false, reason := "h_size_mismatch");
    fi;
    for i in [1..n] do
        if Length(h_mat[i]) <> n then
            return rec(ok := false, reason := "h_size_mismatch");
        fi;
    od;
    if Length(basis_sequence) <> n then
        return rec(ok := false, reason := "basis_size_mismatch");
    fi;
    diag_vals := [];
    for i in [1..n] do
        Add(diag_vals, h_mat[i][i]);
    od;
    alpha_vals := List([1..dim_va], i -> fail);
    beta_vals := List([1..dim_vb], i -> fail);
    for i in [1..n] do
        name := basis_sequence[i];
        if Length(name) >= 2 and name[1] = 'a' then
            idx_str := name{[2..Length(name)]};
            idx := Int(idx_str);
            if idx <> fail and idx >= 1 and idx <= dim_va then
                alpha_vals[idx] := diag_vals[i];
            fi;
        elif Length(name) >= 2 and name[1] = 'b' then
            idx_str := name{[2..Length(name)]};
            idx := Int(idx_str);
            if idx <> fail and idx >= 1 and idx <= dim_vb then
                beta_vals[idx] := diag_vals[i];
            fi;
        fi;
    od;
    good := true;
    for i in [1..dim_va] do
        if alpha_vals[i] = fail then
            good := false;
            break;
        fi;
    od;
    if good then
        for i in [1..dim_vb] do
            if beta_vals[i] = fail then
                good := false;
                break;
            fi;
        od;
    fi;
    delta_f_vals := [];
    if good then
        for i in [1..dim_va] do
            for j in [1..dim_va] do
                if i <> j then
                    Add(delta_f_vals, alpha_vals[i] - alpha_vals[j]);
                fi;
            od;
        od;
        for i in [1..dim_vb] do
            for j in [1..dim_vb] do
                if i <> j then
                    Add(delta_f_vals, beta_vals[i] - beta_vals[j]);
                fi;
            od;
        od;
    fi;
    delta_p_vals := [];
    if good then
        for i in [1..dim_va] do
            for j in [1..dim_vb] do
                Add(delta_p_vals, alpha_vals[i] - beta_vals[j]);
                Add(delta_p_vals, beta_vals[j] - alpha_vals[i]);
            od;
        od;
    fi;
    return rec(
        ok := good,
        diag := diag_vals,
        alpha_values := alpha_vals,
        beta_values := beta_vals,
        delta_f_values := delta_f_vals,
        delta_p_values := delta_p_vals
    );
end;;

BuildSL2RModel := function()
    local Va, Vb, V, dim_va, dim_vb, s_mat, theta, epsilon, omega, bilinear_mat, bilinear_inv, sigma_assoc, sigma_group, sigma_lie, psi_lie, cond_eps, cond_omega, f_desc, p_desc, K_desc, is_in_f, is_in_p, is_in_K, weight_data;
    Va := ["a1"];
    Vb := ["b1"];
    V := Concatenation(Va, Vb);
    dim_va := Length(Va);
    dim_vb := Length(Vb);
    s_mat := [[1, 0], [0, -1]];
    epsilon := 1;
    omega := 1;
    bilinear_mat := [[1, 0], [0, 1]];
    bilinear_inv := bilinear_mat;
    theta := function(x)
        return s_mat * x * s_mat;
    end;
    sigma_assoc := function(x)
        return bilinear_inv * TransposedMat(x) * bilinear_mat;
    end;
    sigma_group := function(g)
        return sigma_assoc(g);
    end;
    sigma_lie := function(x)
        return sigma_assoc(x);
    end;
    psi_lie := function(x)
        return -sigma_assoc(x);
    end;
    is_in_f := function(x)
        return IsInFAIII(x, dim_va, dim_vb);
    end;
    is_in_p := function(x)
        return IsInPAIII(x, dim_va, dim_vb);
    end;
    is_in_K := function(g)
        return is_in_f(g) and DeterminantMat(g) <> 0;
    end;
    cond_eps := (bilinear_mat = epsilon * TransposedMat(bilinear_mat));
    cond_omega := (TransposedMat(s_mat) * bilinear_mat = omega * bilinear_mat * s_mat);
    f_desc := "{X in g | X(V_a)⊆V_a, X(V_b)⊆V_b}";
    p_desc := "{X in g | X(V_a)⊆V_b, X(V_b)⊆V_a}";
    K_desc := "{g in GL(V) | gV_a=V_a, gV_b=V_b}";
    weight_data := BuildAIIIWeightData(dim_va, dim_vb);
    return rec(
        label := "sl(2,R)",
        V := V,
        V_a := Va,
        V_b := Vb,
        epsilon := epsilon,
        omega := omega,
        bilinear_form_matrix := bilinear_mat,
        s := s_mat,
        theta := theta,
        sigma_assoc := sigma_assoc,
        sigma_group := sigma_group,
        sigma_lie := sigma_lie,
        psi_lie := psi_lie,
        checks := rec(epsilon_symmetry := cond_eps, omega_compatibility := cond_omega),
        G_tilde := "GL(V)",
        g_tilde := "Lie(GL(V))=gl(V)",
        G := "{g in GL(V) | sigma(g)^-1=g}",
        g := "{X in gl(V) | psi(X)=-sigma(X)=X}",
        f := f_desc,
        p := p_desc,
        K := K_desc,
        t := weight_data.t_desc,
        alpha_names := weight_data.alpha_names,
        beta_names := weight_data.beta_names,
        alpha_action := weight_data.alpha_action,
        beta_action := weight_data.beta_action,
        Lambda_Va_t := weight_data.Lambda_Va_t,
        Lambda_Vb_t := weight_data.Lambda_Vb_t,
        Delta_g_t := weight_data.Delta_g_t,
        Delta_f_t := weight_data.Delta_f_t,
        Delta_p_t := weight_data.Delta_p_t,
        is_in_f := is_in_f,
        is_in_p := is_in_p,
        is_in_K := is_in_K
    );
end;;

BuildSL2RTriple_Component := function(root_name)
    local model, x, y, h, check_hx, check_hy, check_xy, in_p, in_p_block, basis_sequence, weight_eval;
    model := BuildSL2RModel();
    x := [[0, 1], [0, 0]];
    y := [[0, 0], [1, 0]];
    h := [[1, 0], [0, -1]];
    basis_sequence := ["a1", "b1"];
    weight_eval := EvaluateAIIIWeightsOnH(h, basis_sequence, 1, 1);
    check_hx := MatrixCommutator(h, x) = 2 * x;
    check_hy := MatrixCommutator(h, y) = -2 * y;
    check_xy := MatrixCommutator(x, y) = h;
    in_p := model.theta(x) = -x;
    in_p_block := model.is_in_p(x);
    return rec(
        ok := true,
        root := root_name,
        model := model,
        x := x,
        h := h,
        y := y,
        basis_sequence := basis_sequence,
        weight_eval := weight_eval,
        h_expr := Concatenation("H_{", root_name, "}"),
        e_expr := Concatenation("E_{[", root_name, "]}"),
        f_expr := Concatenation("F_{[", root_name, "]}"),
        checks := rec(hx := check_hx, hy := check_hy, xy := check_xy, x_in_p := in_p, x_in_p_block := in_p_block)
    );
end;;

PrintSL2RTriple_Component := function(root_name)
    local r;
    r := BuildSL2RTriple_Component(root_name);
    if not r.ok then
        Print("        [模块4组件] sl(2,R) 构造跳过\n");
        return r;
    fi;
    Print("        >>> 模块 4(新) 输出: A1 组件构造\n");
    Print("        [模块4组件] V   = ", r.model.V, "\n");
    Print("        [模块4组件] V_a = ", r.model.V_a, "\n");
    Print("        [模块4组件] V_b = ", r.model.V_b, "\n");
    Print("        [模块4组件] ( , ) 矩阵 = ", r.model.bilinear_form_matrix, "\n");
    Print("        [模块4组件] s   = ", r.model.s, "\n");
    Print("        [模块4组件] t   = ", r.model.t, "\n");
    Print("        [模块4组件] alpha = ", r.model.alpha_names, "\n");
    Print("        [模块4组件] beta  = ", r.model.beta_names, "\n");
    Print("        [模块4组件] Delta(g,t) = ", r.model.Delta_g_t, "\n");
    Print("        [模块4组件] Delta(f,t) = ", r.model.Delta_f_t, "\n");
    Print("        [模块4组件] Delta(p,t) = ", r.model.Delta_p_t, "\n");
    if IsRecord(r.weight_eval) and IsBound(r.weight_eval.ok) and r.weight_eval.ok then
        Print("        [模块4组件] h对角(按基顺序) = ", r.weight_eval.diag, "\n");
        Print("        [模块4组件] alpha(h) = ", r.weight_eval.alpha_values, "\n");
        Print("        [模块4组件] beta(h)  = ", r.weight_eval.beta_values, "\n");
        Print("        [模块4组件] Delta(f,t)(h) = ", r.weight_eval.delta_f_values, "\n");
        Print("        [模块4组件] Delta(p,t)(h) = ", r.weight_eval.delta_p_values, "\n");
    fi;
    Print("        [模块4组件] x = ", r.x, "\n");
    Print("        [模块4组件] h = ", r.h, "\n");
    Print("        [模块4组件] y = ", r.y, "\n");
    Print("        [模块4组件] 校验 [h,x]=2x: ", r.checks.hx, "\n");
    Print("        [模块4组件] 校验 [h,y]=-2y: ", r.checks.hy, "\n");
    Print("        [模块4组件] 校验 [x,y]=h: ", r.checks.xy, "\n");
    Print("        [模块4组件] 校验 x∈p(theta(x)=-x): ", r.checks.x_in_p, "\n");
    Print("        [模块4组件] 校验 x∈p(分块定义): ", r.checks.x_in_p_block, "\n");
    Print("        [模块4组件] root = ", root_name, "\n");
    Print("        [模块4组件] h = ", r.h_expr, "\n");
    Print("        [模块4组件] e = ", r.e_expr, "\n");
    Print("        [模块4组件] f = ", r.f_expr, "\n");
    return r;
end;;

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

BuildBDINoncompactCandidates := function(subsystem_roots, colors_opt)
    local candidates, i;
    candidates := [];
    if IsList(colors_opt) and Length(colors_opt) = Length(subsystem_roots) then
        for i in [1..Length(subsystem_roots)] do
            if colors_opt[i] = 1 then
                Add(candidates, subsystem_roots[i]);
            fi;
        od;
    fi;
    if Length(candidates) > 1 then
        Sort(candidates);
    fi;
    return candidates;
end;;

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

SanitizeF4RootList := function(root_list)
    local parser, cleaned, r, v, s;
    cleaned := [];
    parser := fail;
    if IsBound(ParseRootStringFlexible) then
        parser := ParseRootStringFlexible;
    elif IsBound(ParseRootString) then
        parser := ParseRootString;
    fi;
    for r in root_list do
        if parser <> fail then
            v := parser(r);
            if IsList(v) and Length(v) = 4 then
                if IsBoundGlobal("IsF4RootVec") and ValueGlobal("IsF4RootVec")(v) then
                    if IsBoundGlobal("FormatRootF4") then
                        s := ValueGlobal("FormatRootF4")(v);
                    else
                        s := r;
                    fi;
                    Add(cleaned, s);
                fi;
            fi;
        fi;
    od;
    return cleaned;
end;;

CountCharInString := function(s, ch)
    local i, cnt;
    cnt := 0;
    for i in [1..Length(s)] do
        if s[i] = ch then
            cnt := cnt + 1;
        fi;
    od;
    return cnt;
end;;

XiAlpha := function(k, i)
    return Concatenation("alpha_", String(k), "^", String(i));
end;;

XiBeta := function(k, i)
    return Concatenation("beta_", String(k), "^", String(i));
end;;

XiTerm := function(lhs, rhs)
    return Concatenation("X(", lhs, " - ", rhs, ")");
end;;

XiTermsToExpr := function(terms)
    local i, s;
    if Length(terms) = 0 then
        return "0";
    fi;
    s := terms[1];
    for i in [2..Length(terms)] do
        s := Concatenation(s, " + ", terms[i]);
    od;
    return s;
end;;

BuildXiDataForRow := function(row_str, row_idx)
    local n, p, terms, row_type, i, valid_alt, ch_prev, ch_cur;
    n := Length(row_str);
    terms := [];
    if n = 0 then
        return rec(ok := false, reason := "empty_row");
    fi;
    if CountCharInString(row_str, 'a') = n then
        return rec(ok := true, row_type := "all_a", p := 0, terms := [], xi_expr := "0");
    fi;
    if CountCharInString(row_str, 'b') = n then
        return rec(ok := true, row_type := "all_b", p := 0, terms := [], xi_expr := "0");
    fi;
    valid_alt := true;
    ch_prev := row_str[1];
    for i in [2..n] do
        ch_cur := row_str[i];
        if ch_cur = ch_prev then
            valid_alt := false;
            break;
        fi;
        ch_prev := ch_cur;
    od;
    if not valid_alt then
        return rec(ok := false, reason := "row_not_alternating");
    fi;
    if n mod 2 = 0 then
        p := n / 2;
        if row_str[1] = 'b' then
            row_type := "ba...ba";
            for i in Reversed([1..p]) do
                Add(terms, XiTerm(XiAlpha(i, row_idx), XiBeta(i, row_idx)));
                if i > 1 then
                    Add(terms, XiTerm(XiBeta(i, row_idx), XiAlpha(i - 1, row_idx)));
                fi;
            od;
        else
            row_type := "ab...ab";
            for i in Reversed([1..p]) do
                Add(terms, XiTerm(XiBeta(i, row_idx), XiAlpha(i, row_idx)));
                if i > 1 then
                    Add(terms, XiTerm(XiAlpha(i, row_idx), XiBeta(i - 1, row_idx)));
                fi;
            od;
        fi;
    else
        p := (n - 1) / 2;
        if row_str[1] = 'a' then
            row_type := "ab...ba";
            Add(terms, XiTerm(XiAlpha(p + 1, row_idx), XiBeta(p, row_idx)));
            for i in Reversed([1..p]) do
                Add(terms, XiTerm(XiBeta(i, row_idx), XiAlpha(i, row_idx)));
                if i > 1 then
                    Add(terms, XiTerm(XiAlpha(i, row_idx), XiBeta(i - 1, row_idx)));
                fi;
            od;
        else
            row_type := "ba...ab";
            Add(terms, XiTerm(XiBeta(p + 1, row_idx), XiAlpha(p, row_idx)));
            for i in Reversed([1..p]) do
                Add(terms, XiTerm(XiAlpha(i, row_idx), XiBeta(i, row_idx)));
                if i > 1 then
                    Add(terms, XiTerm(XiBeta(i, row_idx), XiAlpha(i - 1, row_idx)));
                fi;
            od;
        fi;
    fi;
    return rec(
        ok := true,
        row_type := row_type,
        p := p,
        terms := terms,
        xi_expr := XiTermsToExpr(terms)
    );
end;;

PrintSL2Triple := function(type, rank, subsystem_roots, colors_opt)
    local engine, ambient;
    engine := "linear_rebuild";
    if IsBoundGlobal("GlobalModule4Engine") then
        engine := ValueGlobal("GlobalModule4Engine");
    fi;
    ambient := "F4";
    if IsBoundGlobal("GlobalAmbientType") then
        ambient := ValueGlobal("GlobalAmbientType");
    fi;
    if ambient = "G2" then
        PrintSL2TripleUnified(type, rank, subsystem_roots, colors_opt);
        return;
    fi;
    PrintSL2TripleUnified(type, rank, subsystem_roots, colors_opt);
end;;

# 解析形如 "sqrt(6)E_{[a1]}" 或 "2H_{a1}" 的项，返回记录 {kind, root, coeff}
ParseTerm := function(term)
    local i, kind, root, coeff_str, coeff, idx_H, idx_E, idx_F, s, start_idx, end_idx, ParseRational, parsed_val;
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

# F4 根向量格式化（简单根基 a1..a4）
FormatRootF4 := function(v)
    local s, i, coeff, term;
    s := "";
    for i in [1..4] do
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

# 判断 F4 非紧根：a1 系数为奇数
IsNoncompactF4Vec := function(v)
    return (v[1] mod 2) <> 0;
end;;

F4_SymMatrix := function()
    return [[4,-2,0,0],[-2,4,-2,0],[0,-2,2,-1],[0,0,-1,2]];
end;;

F4_Cartan := function()
    return [[2,-1,0,0],[-1,2,-2,0],[0,-1,2,-1],[0,0,-1,2]];
end;;

F4_Reflect := function(vec, i)
    local A, j, inner, res;
    A := F4_Cartan();
    inner := 0;
    for j in [1..4] do
        inner := inner + vec[j] * A[j][i];
    od;
    res := ShallowCopy(vec);
    res[i] := res[i] - inner;
    return res;
end;;

F4_GenerateRootSet := function()
    local roots, queue, v, u, exists, i, w;
    roots := [];
    queue := [];
    Add(queue, [1,0,0,0]); Add(queue, [-1,0,0,0]);
    Add(queue, [0,1,0,0]); Add(queue, [0,-1,0,0]);
    Add(queue, [0,0,1,0]); Add(queue, [0,0,-1,0]);
    Add(queue, [0,0,0,1]); Add(queue, [0,0,0,-1]);
    while Length(queue) > 0 do
        v := queue[Length(queue)];
        Unbind(queue[Length(queue)]);
        exists := false;
        for u in roots do
            if u = v then exists := true; break; fi;
        od;
        if not exists then
            Add(roots, v);
            for i in [1..4] do
                w := F4_Reflect(v, i);
                Add(queue, w);
            od;
        fi;
    od;
    return roots;
end;;

IsF4RootVec := function(v)
    local roots, r;
    if not IsList(v) or Length(v) <> 4 then return false; fi;
    if ValueGlobal("GlobalF4RootSet_SimpleCoeffs") = fail then
        GlobalF4RootSet_SimpleCoeffs := ValueGlobal("F4_GenerateRootSet")();
    fi;
    roots := ValueGlobal("GlobalF4RootSet_SimpleCoeffs");
    for r in roots do
        if r = v then return true; fi;
    od;
    return false;
end;;

# 依据颜色列表判断非紧：对白色（非紧）简单根的系数求和取奇偶
IsNoncompactByColorsF4Vec := function(v, colors_list)
    local i, s, c;
    if not IsList(colors_list) or Length(colors_list) < 4 then
        return IsNoncompactF4Vec(v);
    fi;
    s := 0;
    for i in [1..4] do
        if colors_list[i] = 1 then
            c := v[i];
            if c < 0 then c := -c; fi;
            s := s + (c mod 2);
        fi;
    od;
    return (s mod 2) = 1;
end;;
# 构造 BDI 候选非紧正根集：由给定的“非紧基根”产生沿链的连续和
BuildBDINoncompactCandidates := function(subsystem_roots, colors_opt)
    local parser, i, v, name, seen, candidates, preferred, IsNoncompactCandidate, S, Inner, ok, w, j, pv, pw;
    parser := fail;
    if IsBoundGlobal("ParseRootStringFlexible") then
        parser := ValueGlobal("ParseRootStringFlexible");
    elif IsBoundGlobal("ParseRootString") then
        parser := ValueGlobal("ParseRootString");
    fi;
    S := fail;
    if IsBoundGlobal("F4_SymMatrix") then
        S := ValueGlobal("F4_SymMatrix")();
    fi;
    Inner := function(a, b)
        if S = fail then return 0; fi;
        return a[1]*(S[1][1]*b[1]+S[1][2]*b[2]+S[1][3]*b[3]+S[1][4]*b[4])
             + a[2]*(S[2][1]*b[1]+S[2][2]*b[2]+S[2][3]*b[3]+S[2][4]*b[4])
             + a[3]*(S[3][1]*b[1]+S[3][2]*b[2]+S[3][3]*b[3]+S[3][4]*b[4])
             + a[4]*(S[4][1]*b[1]+S[4][2]*b[2]+S[4][3]*b[3]+S[4][4]*b[4]);
    end;
    seen := [];
    candidates := [];
    IsNoncompactCandidate := function(vec)
        if IsBoundGlobal("IsNoncompactF4Vec") then
            return ValueGlobal("IsNoncompactF4Vec")(vec);
        fi;
        return fail;
    end;
    # 预置一条可用的强正交非紧根链，作为后备
    preferred := [
        "a1",
        "-a1-a2-a3",
        "-a1-a2-a3-a4",
        "-a1-a2-2a3-a4",
        "-a1-2a2-2a3-a4",
        "-a1-2a2-3a3-a4"
    ];
    if parser <> fail then
        # 优先使用 subsystem_roots 中的合法非紧根（保持顺序，去重）
        for i in [1..Length(subsystem_roots)] do
            v := parser(subsystem_roots[i]);
            if IsList(v) and Length(v) = 4 then
                if IsF4RootVec(v) and IsNoncompactCandidate(v) then
                    name := FormatRootF4(v);
                    if not name in seen then
                        if name[1] = '-' then
                            if name{[2..Length(name)]} in seen then
                                # skip
                            else
                                ok := true;
                                for j in [1..Length(candidates)] do
                                    pv := parser(candidates[j]); pw := v;
                                    if Inner(pv, pw) <> 0 then ok := false; break; fi;
                                od;
                                if ok then
                                    Add(seen, name);
                                    Add(candidates, name);
                                fi;
                            fi;
                        else
                            if Concatenation("-", name) in seen then
                                # skip
                            else
                                ok := true;
                                for j in [1..Length(candidates)] do
                                    pv := parser(candidates[j]); pw := v;
                                    if Inner(pv, pw) <> 0 then ok := false; break; fi;
                                od;
                                if ok then
                                    Add(seen, name);
                                    Add(candidates, name);
                                fi;
                            fi;
                        fi;
                    fi;
                fi;
            fi;
        od;
        # 其次补充预置强正交链，避免候选不足
        for i in [1..Length(preferred)] do
            v := parser(preferred[i]);
            if IsList(v) and Length(v) = 4 then
                if IsF4RootVec(v) and IsNoncompactCandidate(v) then
                    name := FormatRootF4(v);
                    if not name in seen then
                        if name[1] = '-' then
                            if name{[2..Length(name)]} in seen then
                                # skip
                            else
                                ok := true;
                                for j in [1..Length(candidates)] do
                                    pv := parser(candidates[j]); pw := v;
                                    if Inner(pv, pw) <> 0 then ok := false; break; fi;
                                od;
                                if ok then
                                    Add(seen, name);
                                    Add(candidates, name);
                                fi;
                            fi;
                        else
                            if Concatenation("-", name) in seen then
                                # skip
                            else
                                ok := true;
                                for j in [1..Length(candidates)] do
                                    pv := parser(candidates[j]); pw := v;
                                    if Inner(pv, pw) <> 0 then ok := false; break; fi;
                                od;
                                if ok then
                                    Add(seen, name);
                                    Add(candidates, name);
                                fi;
                            fi;
                        fi;
                    fi;
                fi;
            fi;
        od;
        # 最后从全体 F4 根集中补充更多非紧根（按字典序，去重）
        if (Length(candidates) < 12 and ValueGlobal("GlobalF4RootSet_SimpleCoeffs") <> fail) or IsBoundGlobal("F4_GenerateRootSet") then
            if ValueGlobal("GlobalF4RootSet_SimpleCoeffs") = fail then
                GlobalF4RootSet_SimpleCoeffs := ValueGlobal("F4_GenerateRootSet")();
            fi;
            for v in ValueGlobal("GlobalF4RootSet_SimpleCoeffs") do
                if IsNoncompactCandidate(v) then
                    name := FormatRootF4(v);
                    if not name in seen then
                        ok := true;
                        for j in [1..Length(candidates)] do
                            pv := parser(candidates[j]); pw := v;
                            if Inner(pv, pw) <> 0 then ok := false; break; fi;
                        od;
                        if ok then
                            Add(seen, name);
                            Add(candidates, name);
                        fi;
                    fi;
                fi;
                if Length(candidates) >= 24 then break; fi;
            od;
        fi;
    fi;
    return candidates;
end;;

# 将形如 "sqrt(6)E_{[a1]} + E_{[a1+a2]}" 的字符串解析为 root->系数之和的字典（以列表对实现）
CollectCoeffs := function(expr, expected_kind)
    local parts, p, term, data, acc, found, pair, i, Trim, rest, pos;
    # 简单去除首尾空格
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
VerifyEFH_NormalConditions := function(h_str, e_str, f_str)
    local h_map, e_terms, f_terms, i, n, accum, pair, root, ch, sumv, ok, pr, he_ok, hf_ok, parser, alpha_h, cond, logf, S, h_terms, eval_alpha, has_higher, ParseTermEFH, h_norm, e_norm, f_norm, h_coeffs;
    if IsBoundGlobal("GlobalOutFile") then
        logf := ValueGlobal("GlobalOutFile");
    elif IsBoundGlobal("GlobalOutputDir") then
        logf := Concatenation(ValueGlobal("GlobalOutputDir"), "/F4_pipeline_output.txt");
    else
        logf := "F4_pipeline_output.txt";
    fi;
    h_norm := NormalizeExprForSplitTermsF4(h_str);
    e_norm := NormalizeExprForSplitTermsF4(e_str);
    f_norm := NormalizeExprForSplitTermsF4(f_str);
    # 有些构造会在同一根上出现多次项，此时不能先求和再相乘（会引入交叉项）
    # 采用按序配对的方法：对第 i 项计算 e_i * f_i，并按根名累加
    ParseTermEFH := function(term)
        local s, kind, root, coeff, idx_H, idx_E, idx_F, start_idx, end_idx, i2, coeff_str, eff_coeff, ParseRational, slash, num_str1, num_str2, num1, num2, num_val, parsed_val;
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
        for pr in SplitExprTermsByPlusF4(e_norm) do
            pair := ParseTermEFH(pr);
            if pair.kind = "E" then Add(e_terms, [pair.root, pair.coeff, pair.eff_coeff]); fi;
        od;
    fi;
    if f_norm <> "0" then
        for pr in SplitExprTermsByPlusF4(f_norm) do
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
            AppendTo(logf, Concatenation("        [诊断] 第 ", String(i), " 项根不匹配: E=", e_terms[i][1], " F=", f_terms[i][1], "\n"));
            Print("EFH_CHECK: FAIL\n");
            AppendTo(logf, "EFH_CHECK: FAIL\n");
            return false;
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
        AppendTo(logf, "        校验: [e,f]=h: PASS\n");
        Print("EFH_CHECK: PASS\n");
        AppendTo(logf, "EFH_CHECK: PASS\n");
    else
        Print("        校验: [e,f]=h: FAIL\n");
        AppendTo(logf, "        校验: [e,f]=h: FAIL\n");
        Print("EFH_CHECK: FAIL\n");
        AppendTo(logf, "EFH_CHECK: FAIL\n");
    fi;
    he_ok := true;
    hf_ok := true;
    if (IsBoundGlobal("ParseRootStringFlexible") or IsBoundGlobal("ParseRootString")) then
        has_higher := false;
        if PositionSublist(h_norm, "a3") <> fail or PositionSublist(h_norm, "a4") <> fail then
            has_higher := true;
        elif PositionSublist(e_norm, "a3") <> fail or PositionSublist(e_norm, "a4") <> fail then
            has_higher := true;
        elif PositionSublist(f_norm, "a3") <> fail or PositionSublist(f_norm, "a4") <> fail then
            has_higher := true;
        fi;
        parser := fail;
        if has_higher and IsBoundGlobal("ParseRootStringFlexible") then
            parser := ValueGlobal("ParseRootStringFlexible");
        elif IsBoundGlobal("ParseRootString") then
            parser := ValueGlobal("ParseRootString");
        elif IsBoundGlobal("ParseRootStringFlexible") then
            parser := ValueGlobal("ParseRootStringFlexible");
        fi;
        S := fail;
        if has_higher and IsBoundGlobal("F4_SymMatrix") then
            S := ValueGlobal("F4_SymMatrix")();
        fi;
        h_terms := CollectCoeffs(h_norm, "H");
        h_coeffs := fail;
        if not has_higher and IsBoundGlobal("ParseHString") then
            h_coeffs := ValueGlobal("ParseHString")(h_norm);
        fi;
        eval_alpha := function(alpha_root_str)
            local v, total, pair, g, num, den, j;
            v := parser(alpha_root_str);
            if not has_higher and h_coeffs <> fail and IsList(v) and Length(v) >= 2 then
                if IsBoundGlobal("CalculateWeight") then
                    return ValueGlobal("CalculateWeight")(v, h_coeffs);
                fi;
                return v[1] * (2*h_coeffs[1] - h_coeffs[2]) + v[2] * (-3*h_coeffs[1] + 2*h_coeffs[2]);
            fi;
            total := 0;
            for pair in h_terms do
                g := parser(pair[1]);
                if not has_higher and IsList(v) and IsList(g) and Length(v) >= 2 and Length(g) >= 2 then
                    num := v[1]*(2*g[1]-3*g[2]) + v[2]*(-3*g[1]+6*g[2]);
                    den := g[1]*(2*g[1]-3*g[2]) + g[2]*(-3*g[1]+6*g[2]);
                elif has_higher and IsList(v) and IsList(g) and Length(v) = 4 and Length(g) = 4 and S <> fail then
                    num := v[1]*(S[1][1]*g[1]+S[1][2]*g[2]+S[1][3]*g[3]+S[1][4]*g[4])
                         + v[2]*(S[2][1]*g[1]+S[2][2]*g[2]+S[2][3]*g[3]+S[2][4]*g[4])
                         + v[3]*(S[3][1]*g[1]+S[3][2]*g[2]+S[3][3]*g[3]+S[3][4]*g[4])
                         + v[4]*(S[4][1]*g[1]+S[4][2]*g[2]+S[4][3]*g[3]+S[4][4]*g[4]);
                    den := g[1]*(S[1][1]*g[1]+S[1][2]*g[2]+S[1][3]*g[3]+S[1][4]*g[4])
                         + g[2]*(S[2][1]*g[1]+S[2][2]*g[2]+S[2][3]*g[3]+S[2][4]*g[4])
                         + g[3]*(S[3][1]*g[1]+S[3][2]*g[2]+S[3][3]*g[3]+S[3][4]*g[4])
                         + g[4]*(S[4][1]*g[1]+S[4][2]*g[2]+S[4][3]*g[3]+S[4][4]*g[4]);
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
                AppendTo(logf, Concatenation("        [诊断] E 根 ", e_terms[i][1], " 的 α(h)=", String(alpha_h), " : PASS\n"));
            else
                Print("        [诊断] E 根 ", e_terms[i][1], " 的 α(h)=", alpha_h, " : FAIL\n");
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
                AppendTo(logf, Concatenation("        [诊断] F 根 ", f_terms[i][1], " 的 α(h)=", String(-alpha_h), " => 检查 -α(h)=-2: PASS\n"));
            else
                Print("        [诊断] F 根 ", f_terms[i][1], " 的 α(h)=", -alpha_h, " => 检查 -α(h)=-2: FAIL\n");
                AppendTo(logf, Concatenation("        [诊断] F 根 ", f_terms[i][1], " 的 α(h)=", String(-alpha_h), " => 检查 -α(h)=-2: FAIL\n"));
                hf_ok := false;
                break;
            fi;
        od;
        if he_ok then
            Print("        校验: [h,e]=2e: PASS\n");
            AppendTo(logf, "        校验: [h,e]=2e: PASS\n");
        else
            Print("        校验: [h,e]=2e: FAIL\n");
            AppendTo(logf, "        校验: [h,e]=2e: FAIL\n");
        fi;
        if hf_ok then
            Print("        校验: [h,f]=-2f: PASS\n");
            AppendTo(logf, "        校验: [h,f]=-2f: PASS\n");
        else
            Print("        校验: [h,f]=-2f: FAIL\n");
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

VerifyPK_Subsystem_G2 := function(h_str, e_str, f_str, subsystem_roots, colors_opt)
    local i, root_name, nc_flag, ok, ok_e, ok_f, ok_h, ParseRootSimple, IsNoncompactRoot, CheckTerms, p, data, ComputeSubsystemCoeffs, black_sum;
    ParseRootSimple := ParseRootStringG2Simple;
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
    CheckTerms := function(expr, expect_noncompact)
        local parts, term, nc, NormalizeRootString, root_clean, SplitByPlusOutside, buf, ch, Trim, depth;
        NormalizeRootString := function(s)
            return Filtered(s, c -> c <> ' ' and c <> '[' and c <> ']');
        end;
        Trim := function(s)
            local l, r;
            l := 1; r := Length(s);
            while l <= r and s[l] = ' ' do l := l + 1; od;
            while r >= l and s[r] = ' ' do r := r - 1; od;
            if r < l then return ""; fi;
            return s{[l..r]};
        end;
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

VerifyPK_Subsystem_F4 := function(h_str, e_str, f_str, subsystem_roots, colors_opt)
    local ParseRoot, IsNoncompactRoot, IsNoncompactCandidate, CheckTerms, ok_e, ok_f, pr, pair;
    ParseRoot := fail;
    if IsBoundGlobal("ParseRootStringFlexible") then
        ParseRoot := ValueGlobal("ParseRootStringFlexible");
    elif IsBoundGlobal("ParseRootString") then
        ParseRoot := ValueGlobal("ParseRootString");
    fi;
    IsNoncompactCandidate := function(vec)
        if IsBoundGlobal("IsNoncompactF4Vec") then
            return ValueGlobal("IsNoncompactF4Vec")(vec);
        fi;
        return fail;
    end;
    IsNoncompactRoot := function(root_str)
        local vec;
        if ParseRoot = fail then return fail; fi;
        vec := ParseRoot(root_str);
        if not IsList(vec) or Length(vec) <> 4 then
            return fail;
        fi;
        return IsNoncompactCandidate(vec);
    end;
    CheckTerms := function(expr, expect_noncompact)
        local term, nc;
        if expr = "0" then return true; fi;
        for term in SplitString(expr, " + ") do
            pair := ParseTerm(term);
            if (expect_noncompact and pair.kind = "E") or (expect_noncompact and pair.kind = "F") then
                nc := IsNoncompactRoot(pair.root);
                if nc = fail then
                    Print("        [诊断] P/K 根无法判定: root=", pair.root, "\n");
                    return false;
                fi;
                if nc <> true then
                    Print("        [诊断] P/K 根应非紧但为紧: root=", pair.root, "\n");
                    return false;
                fi;
            fi;
        od;
        return true;
    end;
    ok_e := CheckTerms(e_str, true);
    ok_f := CheckTerms(f_str, true);
    return ok_e and ok_f;
end;;
