class PathologyDiagnosticEngine:
    """
    病理诊断推理引擎
    基于贝叶斯网络 + 决策树 + 规则系统的混合架构
    """
    
    def __init__(self):
        self.clinical_data = {}      # 临床信息
        self.morphology = {}         # 形态学特征
        self.ihc_panel = {}          # 免疫组化结果
        self.molecular = {}          # 分子检测
        self.diagnosis_space = []    # 诊断空间
        self.confidence = 0.0        # 置信度
        
    # ========== 链路1：先验概率计算 ==========
    def calculate_prior(self, clinical_features):
        """
        基于临床特征计算先验概率
        """
        prior_probs = {}
        
        # 年龄因子
        age = clinical_features['age']
        if age > 50:
            prior_probs['DDLPS'] = 0.35
            prior_probs['RMS'] = 0.02  # 成人RMS罕见
        else:
            prior_probs['DDLPS'] = 0.10
            prior_probs['RMS'] = 0.15
            
        # 部位因子
        location = clinical_features['location']
        if location in ['大腿深部', '腹膜后']:
            prior_probs['DDLPS'] *= 2.0  # 好发部位加权
            
        # 尺寸因子
        size = clinical_features['size']
        if size > 15:  # cm
            prior_probs['DDLPS'] *= 1.5  # 巨大肿物倾向
            
        # 归一化
        total = sum(prior_probs.values())
        return {k: v/total for k, v in prior_probs.items()}
    
    # ========== 链路2：形态学模式识别 ==========
    def morphology_analysis(self, histology):
        """
        形态学特征提取与分类
        """
        features = {
            '细胞类型': self.extract_cell_type(histology),
            '异型性': self.grade_atypia(histology),
            '核分裂': self.count_mitosis(histology),
            '坏死': self.detect_necrosis(histology)
        }
        
        # 形态学模式匹配
        if features['细胞类型'] == '梭形细胞':
            if features['异型性'] == '高':
                candidates = ['DDLPS', 'UPS', 'LMS', 'MPNST']
            else:
                candidates = ['WDLPS', '纤维瘤病']
        
        return candidates, features
    
    # ========== 链路3：免疫组化决策树 ==========
    def ihc_decision_tree(self, ihc_results):
        """
        分层免疫组化决策
        """
        # 第一层：谱系归属
        if ihc_results.get('CK') or ihc_results.get('EMA'):
            return self._epithelial_pathway()
        
        # 第二层：间叶分化
        lineage_markers = {
            'smooth_muscle': ['SMA', 'Desmin', 'Caldesmon'],
            'skeletal_muscle': ['Myogenin', 'Desmin', 'MyoD1'],
            'neural': ['S100', 'SOX10'],
            'vascular': ['CD31', 'ERG', 'CD34'],
            'adipocytic': ['MDM2', 'CDK4']
        }
        
        activated_pathways = []
        for lineage, markers in lineage_markers.items():
            if any(ihc_results.get(m, 0) > 0 for m in markers):
                activated_pathways.append(lineage)
        
        # 第三层：特异性诊断
        if 'adipocytic' in activated_pathways:
            return self._liposarcoma_pathway(ihc_results)
        elif 'skeletal_muscle' in activated_pathways:
            return self._rhabdomyosarcoma_pathway(ihc_results)
        else:
            return self._undifferentiated_pathway(ihc_results)
    
    def _liposarcoma_pathway(self, ihc):
        """
        脂肪肉瘤专用通路
        """
        mdm2 = ihc.get('MDM2', 0)
        cdk4 = ihc.get('CDK4', 0)
        
        # 核心逻辑
        if mdm2 > 0 and cdk4 > 0:
            confidence = 0.95
            
            # 检查异源性分化
            if ihc.get('Myoglobin', 0) > 0:
                if ihc.get('Myogenin', 0) == 0:
                    diagnosis = "DDLPS + 异源性肌源性分化"
                else:
                    diagnosis = "DDLPS vs RMS collision tumor"
            else:
                diagnosis = "DDLPS"
            
            return {
                'diagnosis': diagnosis,
                'confidence': confidence,
                'recommendation': 'FISH确认MDM2扩增'
            }
        else:
            return None
    
    def _rhabdomyosarcoma_pathway(self, ihc):
        """
        横纹肌肉瘤专用通路
        """
        myogenin = ihc.get('Myogenin', 0)
        desmin = ihc.get('Desmin', 0)
        
        # 严格标准
        if myogenin > 0 and desmin > 0:
            return {
                'diagnosis': 'RMS (需分型)',
                'confidence': 0.90,
                'subtyping': self._rms_subtype(ihc)
            }
        elif myogenin == 0 and desmin == 0:
            # 强排除逻辑
            return {
                'diagnosis': 'RMS excluded',
                'confidence': 0.98,
                'alternative': '考虑其他梭形细胞肉瘤'
            }
        else:
            return {
                'diagnosis': 'Uncertain',
                'recommendation': '补充MyoD1或分子检测'
            }
    
    # ========== 链路4：贝叶斯网络整合 ==========
    def bayesian_integration(self, prior, likelihood_ihc, likelihood_morphology):
        """
        多证据贝叶斯融合
        """
        posterior = {}
        
        for diagnosis in self.diagnosis_space:
            # 联合似然
            joint_likelihood = (
                likelihood_ihc.get(diagnosis, 0.01) * 
                likelihood_morphology.get(diagnosis, 0.5)
            )
            
            # 后验计算
            posterior[diagnosis] = (
                prior.get(diagnosis, 0.01) * joint_likelihood
            )
        
        # 归一化
        total = sum(posterior.values())
        if total > 0:
            posterior = {k: v/total for k, v in posterior.items()}
        
        return posterior
    
    # ========== 链路5：规则冲突解决 ==========
    def resolve_conflicts(self, evidence_dict):
        """
        当证据冲突时的解决策略
        """
        conflicts = self.detect_conflicts(evidence_dict)
        
        for conflict in conflicts:
            if conflict['type'] == 'MDM2+ but morphology atypical':
                # 分子标志物优先
                resolution = {
                    'priority': 'molecular',
                    'action': '保持DDLPS诊断，注释形态学变异',
                    'confidence_penalty': -0.05
                }
            
            elif conflict['type'] == 'Myoglobin+ but Myogenin-':
                # 逻辑排除优先
                resolution = {
                    'priority': 'exclusion_logic',
                    'interpretation': '异源性分化而非原发RMS',
                    'confidence_boost': +0.10
                }
            
            self.apply_resolution(resolution)
        
        return self.final_diagnosis
    
    # ========== 链路6：不确定性量化 ==========
    def uncertainty_quantification(self):
        """
        诊断不确定性的数学量化
        """
        # 信息熵
        entropy = -sum(
            p * np.log2(p) 
            for p in self.posterior_probs.values() 
            if p > 0
        )
        
        # 置信区间
        top_diagnosis = max(self.posterior_probs, key=self.posterior_probs.get)
        confidence_interval = self.bootstrap_ci(top_diagnosis)
        
        # 敏感性分析
        sensitivity = self.parameter_sensitivity()
        
        return {
            'entropy': entropy,  # 熵越低，诊断越确定
            'CI': confidence_interval,
            'sensitivity': sensitivity,
            'recommendation': self.generate_recommendation(entropy)
        }
    
    def generate_recommendation(self, entropy):
        """
        基于不确定性的后续建议
        """
        if entropy < 0.5:  # 低熵，高确定性
            return "诊断可信，建议临床相关性验证"
        elif 0.5 <= entropy < 1.5:
            return "建议分子检测确认"
        else:  # 高熵，低确定性
            return "建议会诊或多学科讨论"
    
    # ========== 主诊断流程 ==========
    def diagnose(self, case_data):
        """
        完整诊断流程
        """
        # 步骤1：先验概率
        prior = self.calculate_prior(case_data['clinical'])
        
        # 步骤2：形态学分析
        morphology_candidates, morph_features = \
            self.morphology_analysis(case_data['histology'])
        
        # 步骤3：免疫组化决策
        ihc_result = self.ihc_decision_tree(case_data['ihc'])
        
        # 步骤4：贝叶斯整合
        posterior = self.bayesian_integration(
            prior, 
            ihc_result['likelihood'],
            morph_features['likelihood']
        )
        
        # 步骤5：冲突解决
        final = self.resolve_conflicts({
            'prior': prior,
            'posterior': posterior,
            'ihc': ihc_result,
            'morphology': morph_features
        })
        
        # 步骤6：不确定性评估
        uncertainty = self.uncertainty_quantification()
        
        # 步骤7：生成报告
        report = self.generate_diagnostic_report(
            final, uncertainty
        )
        
        return report