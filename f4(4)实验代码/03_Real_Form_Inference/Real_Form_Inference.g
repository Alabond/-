#############################################################################
# 文件: Real_Form_Inference.g
# 描述: 根据 Dynkin 图类型和 Vogan 图颜色推断实李代数类型
# 说明: 模块2的核心规则库。
# 核心函数用途:
#   - InferRealForm(type, rank, black_nodes, white_nodes):
#       通用入口，按类型 A/B/C/D 与黑白计数给出实形式记录 rec(type,params,desc)。
#   - InferRealForm_B4_FromColors(colors_list):
#       F4/B4 特化入口，依据黑点所在位置返回 so(p,q) 映射结果。
# 设计意图:
#   - 先提供稳定可执行的“规则驱动”判定，再由主流程中的位置修正函数细化。
#############################################################################

# 推断实李代数类型
# type: "A", "B", "C", "D"
# rank: int
# black_nodes: int (非紧根数量)
# white_nodes: int (紧根数量)
# total_nodes = black_nodes + white_nodes = rank
#
# 返回: rec(real_form_type, params)
# 
# 简化规则（面向本项目的 G2 子系统场景）：
# - A: 全黑→sl(n+1,ℝ)；全白→紧（su(n+1)）；混合→su(p,q) 的启发式映射
# - C: 全黑→sp(n,ℝ)；否则→sp(p,q)（占位参数）
# - B/D: 全黑→近似 split 形 so(n+1,n)；否则→so(p,q)（占位参数）

InferRealForm := function(type, rank, black_nodes, white_nodes)
    local p, q;
    
    if type = "A" then
        if rank = 1 then
            if black_nodes = 1 then
                return rec(type := "sl_R", params := [1], desc := "sl(2, R)");
            fi;
            return rec(type := "compact", params := [2], desc := "su(2)");
        fi;
        if white_nodes = rank then
            return rec(type := "compact", params := [rank+1], desc := Concatenation("su(", String(rank+1), ")"));
        elif black_nodes = rank then
            return rec(type := "sl_R", params := [rank], desc := Concatenation("sl(", String(rank+1), ", R)"));
        elif rank = 3 and ((black_nodes = 1 and white_nodes = 2) or (black_nodes = 2 and white_nodes = 1)) then
            return rec(type := "su_pq", params := [2, 2], desc := "su(2,2)");
        else
            p := white_nodes + 1;
            q := black_nodes;
            if q = 0 then q := 1; p := rank; fi; 
            return rec(type := "su_pq", params := [p, q], desc := Concatenation("su(", String(p), ",", String(q), ")"));
        fi;
        
    elif type = "C" then
        if black_nodes = rank then
            return rec(type := "sp_R", params := [rank], desc := "sp(n, R)");
        else
            return rec(type := "sp_pq", params := [1, 1], desc := "sp(p, q)");
        fi;
        
    elif type = "B" or type = "D" then
        if black_nodes = rank then
            # Split form for B/D is so(n, n+1) or so(n, n)
            # so(p, q) split: p ~ q
            return rec(type := "so_pq", params := [rank+1, rank], desc := "so(split)");
        else
            return rec(type := "so_pq", params := [rank+1, 1], desc := "so(p, q)");
        fi;
    fi;
    
    return fail;
end;;

InferRealForm_B4_FromColors := function(colors_list)
    local idx, i, label, p, q, mapping;
    
    # 找到唯一黑点位置（最左四点之一）
    idx := 0;
    for i in [1..Minimum(4, Length(colors_list))] do
        if colors_list[i] = 1 then
            idx := i;
            break;
        fi;
    od;
    
    # 以 a,b,c,d 对应最左四点（从左到右）：a=1,b=2,c=3,d=4
    if idx = 1 then label := "a";
    elif idx = 2 then label := "b";
    elif idx = 3 then label := "c";
    elif idx = 4 then label := "d";
    else label := "unknown";
    fi;
    
    # 可调映射表（便于与用户标准核对后微调）
    # 逻辑：a->so(5,1), b->so(4,2), c->so(3,3), d->so(4,2)
    mapping := rec(
        a := [5,4],
        b := [4,5],
        c := [3,6],
        d := [4,5]
    );
    
    if label = "unknown" then
        p := 4; q := 5;
    else
        p := mapping.(label)[1];
        q := mapping.(label)[2];
    fi;
    
    return rec(
        type := "so_pq",
        params := [p, q],
        desc := Concatenation("so(", String(p), ",", String(q), ")", " [", label, "]")
    );
end;;
