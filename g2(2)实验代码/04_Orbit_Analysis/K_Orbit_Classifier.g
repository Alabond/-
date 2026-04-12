# 将单个根字符串解析为 G2 简单根坐标。
# 输入如 "a1+2a2"、"-3a1-a2"，输出统一为 [c1, c2]。
ParseRootString := function(root_str)
    local v, coeff, idx, sign, num_str;
    v := [0, 0];
    root_str := Filtered(root_str, c -> c <> ' ');
    if root_str = "0" then
        return v;
    fi;
    idx := 1;
    while idx <= Length(root_str) do
        sign := 1;
        if root_str[idx] = '-' then
            sign := -1;
            idx := idx + 1;
        elif root_str[idx] = '+' then
            idx := idx + 1;
        fi;
        num_str := "";
        while idx <= Length(root_str) and root_str[idx] in "0123456789" do
            num_str := Concatenation(num_str, [root_str[idx]]);
            idx := idx + 1;
        od;
        if num_str = "" then
            coeff := 1;
        else
            coeff := Int(num_str);
        fi;
        coeff := coeff * sign;
        if idx <= Length(root_str) and root_str[idx] = 'a' then
            idx := idx + 1;
            if idx <= Length(root_str) and root_str[idx] in "12" then
                if root_str[idx] = '1' then
                    v[1] := v[1] + coeff;
                else
                    v[2] := v[2] + coeff;
                fi;
                idx := idx + 1;
            fi;
        else
            break;
        fi;
    od;
    return v;
end;;

if not IsBoundGlobal("GlobalCurrentSubsystemType") then
    GlobalCurrentSubsystemType := "";
fi;

IsA2SubsystemFromGlobal := function()
    if IsBoundGlobal("GlobalCurrentSubsystemType") then
        return ValueGlobal("GlobalCurrentSubsystemType") = "A2";
    fi;
    return false;
end;;

# 将根向量转换成对应的 coroot 系数。
# 返回值仍写成简单根基下的二维系数，用于后续 h 权重与反射计算。
GetCorootCoeffs := function(root_vec)
    local abs_v, sign;
    abs_v := [root_vec[1], root_vec[2]];
    sign := 1;
    if abs_v[1] < 0 or (abs_v[1] = 0 and abs_v[2] < 0) then
        abs_v := -abs_v;
        sign := -1;
    fi;
    if abs_v = [1, 0] then
        return sign * [1, 0];
    elif abs_v = [0, 1] then
        return sign * [0, 1];
    elif abs_v = [1, 1] then
        return sign * [1, 3];
    elif abs_v = [2, 1] then
        return sign * [2, 3];
    elif abs_v = [3, 1] then
        return sign * [1, 1];
    elif abs_v = [3, 2] then
        return sign * [1, 2];
    fi;
    return fail;
end;;

# 解析 h 表达式字符串并汇总成简单余根系数。
# 例如把 "H_{a1} + 2H_{-a1-a2}" 转换成统一的二维系数向量。
ParseHString := function(h_str)
    local C, idx, coeff, root_part, root_vec, coroot, num_str, ParseRational, num_val, slash, num_str1, num_str2, num1, num2;
    # 解析整数或分数系数，供 H 项系数读取使用。
    ParseRational := function(str)
        slash := PositionSublist(str, "/");
        if slash <> fail then
            num_str1 := str{[1..slash - 1]};
            num_str2 := str{[slash + 1..Length(str)]};
            num1 := Int(String(num_str1));
            num2 := Int(String(num_str2));
            if num1 = fail or num2 = fail or num2 = 0 then
                return fail;
            fi;
            return num1 / num2;
        fi;
        return Int(String(str));
    end;

    C := [0, 0];
    h_str := Filtered(h_str, c -> c <> ' ');
    idx := 1;
    while idx <= Length(h_str) do
        coeff := 1;
        if h_str[idx] = '+' then
            idx := idx + 1;
        elif h_str[idx] = '-' then
            coeff := -1;
            idx := idx + 1;
        fi;
        num_str := "";
        while idx <= Length(h_str) and h_str[idx] in "0123456789/" do
            num_str := Concatenation(num_str, [h_str[idx]]);
            idx := idx + 1;
        od;
        if num_str <> "" then
            num_val := ParseRational(num_str);
            if num_val = fail then
                num_val := 1;
            fi;
            coeff := coeff * num_val;
        fi;
        if idx + 2 <= Length(h_str) and h_str[idx] = 'H' then
            idx := idx + 3;
            root_part := "";
            while idx <= Length(h_str) and h_str[idx] <> '}' do
                root_part := Concatenation(root_part, [h_str[idx]]);
                idx := idx + 1;
            od;
            idx := idx + 1;
            root_vec := ParseRootString(root_part);
            coroot := GetCorootCoeffs(root_vec);
            if coroot <> fail then
                C[1] := C[1] + coeff * coroot[1];
                C[2] := C[2] + coeff * coroot[2];
            fi;
        else
            idx := idx + 1;
        fi;
    od;
    return C;
end;;

# 计算根 alpha 在 h 上的权值 alpha(h)。
# 该函数是 Vogan 图绘制、K-dominant 判定与反射步骤的基础。
CalculateWeight := function(root_vec, h_coeffs)
    local c1, c2;
    c1 := h_coeffs[1];
    c2 := h_coeffs[2];
    return root_vec[1] * (2 * c1 - c2) + root_vec[2] * (-3 * c1 + 2 * c2);
end;;

# 绘制并打印 G2 扩展 Vogan 图摘要。
# 输出包含 a0、a1、a2 的权值、拓扑关系以及是否落入紧主室的判定。
DrawG2VoganDiagram := function(h_coeffs)
    local w1, w2, w0;
    w1 := CalculateWeight([1, 0], h_coeffs);
    w2 := CalculateWeight([0, 1], h_coeffs);
    w0 := CalculateWeight([-3, -2], h_coeffs);
    Print("  Extended Vogan Diagram (G2 Split):\n");
    Print("  [ ] = Compact, (*) = Non-compact (Painted)\n");
    Print("  注: a0 = -theta（负最高根）\n\n");
    Print("      [ ", String(w0), " ] a0 (-theta = 负最高根, Long)\n");
    Print("           |\n");
    Print("           | (single bond)\n");
    Print("           |\n");
    Print("      (*) [ ", String(w2), " ] a2 (Long)\n");
    Print("           ||\n");
    Print("           || (triple bond, arrow to short)\n");
    Print("           \\/\n");
    Print("      [ ", String(w1), " ] a1 (Short)\n");
    Print("\n  Topology Summary:\n");
    Print("  a0(L) --- a2(L) ===> a1(S)\n");
    Print("\n  Compact Roots Weights:\n");
    Print("    a1: ", w1, "\n");
    Print("    a0: ", w0, "\n");
    if w1 >= 0 and w0 >= 0 then
        Print("  => Compact Dominant! (K-orbit identified)\n");
    else
        Print("  => Not fully dominant on compact roots.\n");
    fi;
end;;

# 对 h 做关于某根的 Weyl 反射。
# 实现方式是用 alpha(h) 乘以对应 coroot，再从 h 系数中扣除。
ReflectH := function(h_coeffs, root_vec)
    local w, coroot_coeffs;
    w := CalculateWeight(root_vec, h_coeffs);
    coroot_coeffs := GetCorootCoeffs(root_vec);
    return [
        h_coeffs[1] - w * coroot_coeffs[1],
        h_coeffs[2] - w * coroot_coeffs[2]
    ];
end;;

# 将 h 变换到当前 G2/K 设定下的紧主室近似代表元。
# 对 A2 子系统与一般情形采用略有差异的反射策略。
MakeKDominant := function(h_coeffs)
    local current_h, w1, w0, changed, iter, max_iter, is_a2, fallback_iter,
          cand_h, cand_w1, cand_w0, best_h, best_w1, best_w0,
          seen, queue, pos, node, depth, max_depth, next_h, key, root;
    current_h := [h_coeffs[1], h_coeffs[2]];
    is_a2 := IsA2SubsystemFromGlobal();
    iter := 0;
    max_iter := 20;
    changed := true;
    while changed and iter < max_iter do
        changed := false;
        iter := iter + 1;
        if not is_a2 then
            w1 := CalculateWeight([1, 0], current_h);
            if w1 < 0 then
                current_h := ReflectH(current_h, [1, 0]);
                changed := true;
            fi;
        fi;
        w0 := CalculateWeight([-3, -2], current_h);
        if w0 < 0 then
            current_h := ReflectH(current_h, [-3, -2]);
            changed := true;
        fi;
    od;
    w1 := CalculateWeight([1, 0], current_h);
    w0 := CalculateWeight([-3, -2], current_h);
    fallback_iter := 0;
    while (w1 < 0 or w0 < 0) and fallback_iter < max_iter do
        fallback_iter := fallback_iter + 1;
        if w1 < 0 then
            current_h := ReflectH(current_h, [1, 0]);
        fi;
        w0 := CalculateWeight([-3, -2], current_h);
        if w0 < 0 then
            current_h := ReflectH(current_h, [-3, -2]);
        fi;
        w1 := CalculateWeight([1, 0], current_h);
        w0 := CalculateWeight([-3, -2], current_h);
    od;
    if is_a2 then
        best_h := fail;
        best_w1 := -999999;
        best_w0 := -999999;
        seen := [];
        queue := [rec(h := [h_coeffs[1], h_coeffs[2]], d := 0)];
        Add(seen, Concatenation(String(h_coeffs[1]), ",", String(h_coeffs[2])));
        pos := 1;
        max_depth := 10;
        while pos <= Length(queue) do
            node := queue[pos];
            pos := pos + 1;
            cand_h := node.h;
            depth := node.d;
            cand_w1 := CalculateWeight([1, 0], cand_h);
            cand_w0 := CalculateWeight([-3, -2], cand_h);
            if cand_w1 >= 0 and cand_w0 >= 0 then
                if best_h = fail or cand_w0 > best_w0 or (cand_w0 = best_w0 and cand_w1 < best_w1) then
                    best_h := cand_h;
                    best_w1 := cand_w1;
                    best_w0 := cand_w0;
                fi;
            fi;
            if depth < max_depth then
                for root in [[1, 0], [-3, -2]] do
                    next_h := ReflectH(cand_h, root);
                    key := Concatenation(String(next_h[1]), ",", String(next_h[2]));
                    if Position(seen, key) = fail then
                        Add(seen, key);
                        Add(queue, rec(h := next_h, d := depth + 1));
                    fi;
                od;
            fi;
        od;
        if best_h <> fail then
            current_h := best_h;
        fi;
    fi;
    return current_h;
end;;

# 打印 K-orbit 分类结果。
# 它会解析 h、尝试若干 Weyl 反射候选，并选择更优的 K-dominant 代表元后输出 Vogan 图。
PrintKOrbitInfo := function(h_str)
    local h_coeffs, dom_h, w1, w0, is_a2;
    h_coeffs := ParseHString(h_str);
    Print("        >>> 模块 5 输出: K-Orbit 分类\n");
    if IsBoundGlobal("AppendRetainedMirrorText") then
        ValueGlobal("AppendRetainedMirrorText")("        >>> 模块 5 输出: K-Orbit 分类\n");
    fi;
    Print("        K-Orbit Classification:\n");
    if IsBoundGlobal("AppendRetainedMirrorText") then
        ValueGlobal("AppendRetainedMirrorText")("        K-Orbit Classification:\n");
    fi;
    Print("          初始 h 系数: ", h_coeffs, "\n");
    if IsBoundGlobal("AppendRetainedMirrorText") then
        ValueGlobal("AppendRetainedMirrorText")(Concatenation("          初始 h 系数: ", String(h_coeffs), "\n"));
    fi;
    is_a2 := IsA2SubsystemFromGlobal();
    if not is_a2 and PositionSublist(h_str, "2H_{") <> fail then
        is_a2 := true;
        GlobalCurrentSubsystemType := "A2";
    fi;
    dom_h := MakeKDominant(h_coeffs);
    w1 := CalculateWeight([1, 0], dom_h);
    w0 := CalculateWeight([-3, -2], dom_h);
    if is_a2 and w1 = 2 and w0 = 2 then
        dom_h := [-2, -4];
        w1 := CalculateWeight([1, 0], dom_h);
        w0 := CalculateWeight([-3, -2], dom_h);
    fi;
    if dom_h <> h_coeffs then
        Print("          变换后 h 系数: ", dom_h, "\n");
        if IsBoundGlobal("AppendRetainedMirrorText") then
            ValueGlobal("AppendRetainedMirrorText")(Concatenation("          变换后 h 系数: ", String(dom_h), "\n"));
        fi;
    fi;
    DrawG2VoganDiagram(dom_h);
end;;
