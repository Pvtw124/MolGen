def prepare_rigid_crossover(_smiles, side, consider_ring=True, _minimum_len=4):
    """ 1 point crossover
    :param _smiles: SMILES (str)
    :param side: Left SMILES or Right SMILES ['L'|'R'] (str)
    :param consider_ring: consider about avoiding ring (bool)
    :param _minimum_len: minimum cut size (int)
    :return:
    """

    if side not in ["L", "R"]:
        raise Exception("You must choice in L(Left) or R(Right)")

    _smiles_len = len(_smiles)
    _smi = None

    if consider_ring:  # ring 부분을 피해서 자르도록 고려.
        avoid_ring_list = set_avoid_ring(_smiles)

    p = 0
    _start = None
    _end = None
    _gate = False
    while not _gate:  # 원하는 형태의 smiles piece가 나올 때 까지 반복
        if p == limit_:
            raise ValueError(f"main_gate fail ({side}): {_smiles}")

        if consider_ring:
            j = 0
            ring_gate = False
            if side == "L":
                while not ring_gate:
                    if j == 30:
                        raise ValueError(f"ring_gate fail (L): {_smiles}")
                    _end = np.random.randint(_minimum_len, _smiles_len + 1)
                    if _end not in avoid_ring_list:
                        ring_gate = True
                    j += 1
            elif side == "R":
                while not ring_gate:
                    if j == 30:
                        raise ValueError(f"ring_gate fail (R): {_smiles}")
                    _start = np.random.randint(0, _smiles_len - _minimum_len)
                    if _start not in avoid_ring_list:
                        ring_gate = True
                    j += 1
            _smi = _smiles[_start:_end]
        else:
            if side == "L":
                _end = np.random.randint(_minimum_len, _smiles_len)
            elif side == "R":
                _start = np.random.randint(0, _smiles_len - _minimum_len)

            _smi = _smiles[_start:_end]
            chk_ring = re.findall(r"\d", _smi)
            i = 0
            for i in set(chk_ring):
                list_ring = [_ for _, val in enumerate(_smi) if val == i]
                if (len(list_ring) % 2) == 1:
                    b = random.sample(list_ring, 1)
                    _smi = _smi[:b[0]] + _smi[b[0] + 1:]
                    # print(f'@ {_smi} // {_smiles}')

        p += 1

        if "." in _smi:  # 이온은 패스한다.
            continue

        n_chk = 0
        for j in _smi:  # [] 닫혀 있는지 확인.
            if j == "[":
                n_chk += 1
            if j == "]":
                n_chk -= 1
        if n_chk == 0:
            _gate = True

    return _smi
def cut_smi(smi1, smi2, func, ring_bool):

    l_smi = None
    r_smi = None

    try:
        l_smi = func(smi1, "L", ring_bool, 4)
        r_smi = func(smi2, "R", ring_bool, 4)
    except (IndexError, ValueError):
        fail_f.write(f"{l_smi},{r_smi},piece\n")
        raise PermissionError

    return l_smi, r_smi



def get_sliced_smiles(smi1, smi2, func, ring_bool):
    l_smi = None
    r_smi = None

    gate = 0
    while not (l_smi and r_smi):
        gate += 1
        if gate > 10:
            # fail_f.write(f"{l_smi},{r_smi},np\n")
            raise PermissionError
        try:
            l_smi, r_smi = cut_smi(smi1, smi2, func, ring_bool)
        except:
            pass
    return l_smi, r_smi
    
def tight_rm_branch(_smi_l, _smi_r):
    # tmp = time.time()

    _new_smi = _smi_l + _smi_r

    open_branch = get_open_branch(_new_smi)
    close_branch = get_close_branch(_new_smi)

    b = None
    n_branch = chk_branch(_new_smi)

    q = len(_smi_l)
    while n_branch > 0:  # over opened-branch
        _smi_l_open_branch = get_open_branch(_smi_l)
        _smi_r_open_branch = get_open_branch(_smi_r)
        open_branch = get_open_branch(_smi_l + _smi_r)
        avoid_tokens = [
            i for i, e in enumerate(_smi_l + _smi_r)
            if e in ["=", "#", "@", "1", "2", "3", "4", "5", "6", "7", "8"]
        ]

        if len(_smi_r_open_branch) == 0:  # open branch 가 없을 경우
            _smi_r_open_branch.append(len(_smi_r))
        if len(_smi_l_open_branch) == 0:
            _smi_l_open_branch.append(0)

        n = np.random.rand()  # 임의적으로 close branch 를 추가하거나 제거한다.
        if n > 0.5:  # 추가
            branch_gate = False
            j = 0
            while not branch_gate:  # Ring 부분을 피해서 자름
                if j == limit_:
                    raise ValueError
                b = np.random.randint(_smi_l_open_branch[-1] + 1,
                                      _smi_r_open_branch[-1] + q)
                j += 1
                if b not in avoid_tokens:
                    branch_gate = True
            n_branch -= 1
            if b <= len(_smi_l
                        ):  # SMILES 길이를 고려하여 자른다. 좌측 SMILES의 open branch를 cut!
                _smi_l = _smi_l[:b] + ")" + _smi_l[b:]
                q += 1
            else:  # 좌측 SMILES 길이를 제외한 수가 우측 SMILES 문자의 위치를 의미한다.
                b -= len(_smi_l)
                _smi_r = _smi_r[:b] + ")" + _smi_r[b:]
        else:  # 제거
            b = _smi_l_open_branch[-1]  # (Random으로도 가능함. 과한 부분만 Cut!)
            n_branch -= 1
            q -= 1
            _smi_l = _smi_l[:b] + _smi_l[b + 1:]

    while n_branch < 0:  # over closed-branch
        _smi_l_close_branch = get_close_branch(_smi_l)
        _smi_r_close_branch = get_close_branch(_smi_r)
        close_branch = get_close_branch(_smi_l + _smi_r)
        avoid_tokens = [
            i for i, e in enumerate(_smi_l + _smi_r)
            if e in ["=", "#", "@", "1", "2", "3", "4", "5", "6", "7", "8"]
        ]

        if len(_smi_r_close_branch) == 0:
            _smi_r_close_branch.append(len(_smi_r))
        if len(_smi_l_close_branch) == 0:
            _smi_l_close_branch.append(0)

        n = np.random.rand()
        if n > 0.5:
            branch_gate = False
            j = 0
            while not branch_gate:  # Ring 부분을 피해서 자름
                b = np.random.randint(_smi_l_close_branch[-1] + 1,
                                      _smi_r_close_branch[0] + q + 1)
                j += 1
                if b not in (close_branch + avoid_tokens):
                    branch_gate = True
                if j == limit_:
                    raise ValueError
            n_branch += 1
            if b < len(_smi_l):
                _smi_l = _smi_l[:b] + "(" + _smi_l[b:]
                q += 1
            else:
                b -= len(_smi_l)
                _smi_r = _smi_r[:b] + "(" + _smi_r[b:]
        else:
            b = _smi_r_close_branch[0]
            n_branch += 1
            # print(f'{_smi_r[b]}')
            _smi_r = _smi_r[:b] + _smi_r[b + 1:]

    # time_.append(time.time() - tmp)

    return _smi_l + _smi_r



def crossover_smiles(smi1, smi2, func, ring_bool):
    new_smi = None
    mol = None

    l_smi, r_smi = get_sliced_smiles(smi1, smi2, func, ring_bool)

    gate = 0
    while not mol:
        gate += 1
        if gate > 5:
            break
        try:
            new_smi = ModSMI.tight_rm_branch(l_smi, r_smi)
        except ValueError:
            continue
        mol = Chem.MolFromSmiles(new_smi)

    if not mol:
        l_smi, r_smi = get_sliced_smiles(smi2, smi1, func, ring_bool)

        gate = 0
        while not mol:
            gate += 1
            if gate > 5:
                break
            try:
                new_smi = ModSMI.tight_rm_branch(l_smi, r_smi)
            except ValueError:
                continue
            mol = Chem.MolFromSmiles(new_smi)

    return new_smi, mol



new_smi, mol = crossover_smiles(smi1, smi2, ModSMI.prepare_rigid_crossover, True)