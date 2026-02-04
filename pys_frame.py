# 定义诊断空间
Ω = {DDLPS, RMS, UPS, LMS, 其他肉瘤}

# 观察证据集
E = {
    临床: (年龄=59, 性别=女, 部位=大腿深部, 大小=20cm),
    形态: (梭形细胞, 高异型性, 核分裂多, 多核巨细胞),
    IHC: {
        MDM2: 1,      # 阳性=1
        CDK4: 0.3,    # 少量阳性=0.3
        Myogenin: 0,  # 阴性=0
        Desmin: 0,
        Myoglobin: 1,
        CD99: 1,
        其他: 全阴性
    }
}

# 先验概率（基于临床特征）
P(DDLPS | 临床) = 0.35  # 大腿深部巨大肿物，中老年
P(RMS | 临床) = 0.05    # 成人多形性RMS罕见
P(UPS | 临床) = 0.30    # 常见高级别肉瘤
P(LMS | 临床) = 0.20    # 平滑肌肉瘤
P(其他 | 临床) = 0.10

# 似然函数（证据在各诊断下的概率）
P(E_IHC | DDLPS) = P(MDM2=1|DDLPS) × P(CDK4>0|DDLPS) × 
                   P(Myogenin=0|DDLPS) × P(Myoglobin=1|DDLPS异源性)
                 = 0.97 × 0.92 × 0.95 × 0.08
                 ≈ 0.068

P(E_IHC | RMS) = P(MDM2=1|RMS) × P(Myogenin=0|RMS)
               ≈ 0.01 × 0.05  # MDM2在RMS中罕见，Myogenin阴性罕见
               ≈ 0.0005

P(E_IHC | UPS) = P(MDM2=1|UPS) × ...
               ≈ 0.02  # UPS中MDM2罕见阳性

# 贝叶斯更新
P(DDLPS | E) = P(E_IHC | DDLPS) × P(DDLPS | 临床) / P(E_IHC)
             = (0.068 × 0.35) / Σ[P(E_IHC|dᵢ) × P(dᵢ|临床)]
             = 0.0238 / (0.0238 + 0.000025 + 0.006 + ...)
             ≈ 0.0238 / 0.032
             ≈ 0.74

# 加入形态学证据后
P(DDLPS | E_完整) ≈ 0.95
```

### **2.2 决策树数学模型**
```
定义决策函数 D: E → Ω

D(E) = argmax_{d∈Ω} [P(d|E) × U(d)]

其中 U(d) 为效用函数（诊断正确的临床价值）

决策规则：
IF P(DDLPS|E) > θ₁ = 0.90 THEN 诊断=DDLPS
ELSE IF P(DDLPS|E) > θ₂ = 0.70 THEN 建议分子检测
ELSE 重新评估

实际计算：
P(DDLPS|E_本例) ≈ 0.95 > θ₁
→ 诊断成立，但建议分子确认（穿刺样本限制）
```

### **2.3 逻辑代数表达**
```
定义命题变量：
M := MDM2阳性
C := CDK4阳性  
G := Myogenin阳性
D := Desmin阳性
Y := Myoglobin阳性
R := 诊断为RMS
L := 诊断为DDLPS

逻辑公式：

1. DDLPS充分条件：
   (M ∧ C) → L  [置信度 0.95]

2. RMS必要条件：
   R → (G ∨ D)  [逆否命题]
   (¬G ∧ ¬D) → ¬R  [强排除]

3. 异源性分化解释：
   (L ∧ Y) → "DDLPS伴异源性肌源性分化"
   
4. 完整诊断逻辑：
   [(M ∧ C) ∧ (¬G ∧ ¬D) ∧ Y] → "DDLPS + 异源性分化"

真值赋值（本例）：
M=T, C=T, G=F, D=F, Y=T
→ L=T (DDLPS成立)
→ R=F (RMS排除)