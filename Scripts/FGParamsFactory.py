import copy

class FGParamsFactory:
  def __init__(self,
               res_from= None,
               res_to= None,
               self_k= None,
               self_range= None,
               kap_k= None,
               kap_range= None,
               nonspec_k= None,
               nonspec_range= None,
               backbone_k= None,
               backbone_tau= None):
#    assert(res_to>=res_from)
    self.res_from= res_from
    self.res_to= res_to
    self.self_k= self_k
    self.self_range= self_range
    self.kap_k= kap_k
    self.kap_range= kap_range
    self.nonspec_k= nonspec_k
    self.nonspec_range= nonspec_range
    self.backbone_k= backbone_k
    self.backbone_tau= backbone_tau

  def get_copy(self, res_from=None, res_to=None, self_k=None, nonspec_k=None):
    ret= copy.deepcopy(self)
    ret.res_from= res_from
    ret.res_to= res_to
    if self_k is not None:
      ret.self_k= self_k
    if nonspec_k is not None:
      ret.nonspec_k= nonspec_k
    return ret
