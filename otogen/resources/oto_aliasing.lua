local hiragana_map_str = [[
息,bre
あ,a
い,i
う,u
え,e
お,o
ん,n
か,ka
き,ki
く,ku
け,ke
こ,ko
さ,sa
すぃ,si
す,su
せ,se
そ,so
しゃ,sha
し,shi
しゅ,shu
しぇ,she
しょ,sho
ちゃ,cha
ち,chi
ちゅ,chu
ちぇ,che
ちょ,cho
つぁ,tsa
つぃ,tsi
つ,tsu
つぇ,tse
つょ,tso
た,ta
てぃ,ti
てぅ,tu
て,te
と,to
な,na
に,ni
ぬ,nu
ね,ne
の,no
は,ha
ひ,hi
ふ,hu
へ,he
ほ,ho
ふぁ,fa
ふぃ,fi
ふぅ,fu
ふぇ,fe
ふぉ,fo
ま,ma
み,mi
む,mu
め,me
も,mo
や,ya
やぃ,yi
ゆ,yu
やぇ,ye
よ,yo
ら,ra
り,ri
る,ru
れ,re
ろ,ro
らぁ,la
りぃ,li
るぅ,lu
れぇ,le
ろぉ,lo
らら,rra
りり,rri
るる,rru
れれ,rre
ろろ,rro
わ,wa
を,wo
うぃ,wi
うぅ,wu
うぇ,we
うぉ,wo
が,ga
ぎ,gi
ぐ,gu
げ,ge
ご,go
じぁ,ja
じ,ji
じぅ,ju
ぢぇ,je
じぉ,jo
じゃ,jya
じぃ,jyi
じゅ,jyu
じぇ,jye
じょ,jyo
ざ,za
ずぃ,zi
ず,zu
ぜ,ze
ぞ,zo
だ,da
でぃ,di
でぅ,du
で,de
ど,do
ば,ba
び,bi
ぶ,bu
べ,be
ぼ,bo
ヴぁ,va
ヴぃ,vi
ヴ,vu
ヴぇ,ve
ヴぉ,vo
ぱ,pa
ぴ,pi
ぷ,pu
ぺ,pe
ぽ,po
きゃ,kya
きぃ,kyi
きゅ,kyu
きぇ,kye
きょ,kyo
てゃ,tya
てぃ,tyi
てゅ,tyu
てぃぇ,tye
てょ,tyo
にゃ,nya
にぃ,nyi
にゅ,nyu
にぇ,nye
にょ,nyo
ひゃ,hya
ひぃ,hyi
ひゅ,hyu
ひぇ,hye
ひょ,hyo
みゃ,mya
みぃ,myi
みゅ,myu
みぇ,mye
みょ,myo
りゃ,rya
りぃ,ryi
りゅ,ryu
りぇ,rye
りょ,ryo
りりゃ,rrya
りりぃ,rryi
りりゅ,rryu
りりぇ,rrye
りりょ,rryo
ぎゃ,gya
ぎぃ,gyi
ぎゅ,gyu
ぎぇ,gye
ぎょ,gyo
びゃ,bya
びぃ,byi
びゅ,byu
びぇ,bye
びょ,byo
ぴゃ,pya
ぴぃ,pyi
ぴゅ,pyu
ぴぇ,pye
ぴょ,pyo
]]
local phoneme_map_str = [[
rryo,ry,o
a,a
i,i
u,u
e,e
o,o
n,N
ka,k,a
ki,k,i
ku,k,u
ke,k,e
ko,k,o
sa,s,a
si,s,i
su,s,u
se,s,e
so,s,o
sha,sh,a
shi,sh,i
shu,sh,u
she,sh,e
sho,sh,o
cha,ch,a
chi,ch,i
chu,ch,u
che,ch,e
cho,ch,o
tsa,ts,a
tsi,ts,i
tsu,ts,u
tse,ts,e
tso,ts,o
ta,t,a
ti,t,i
tu,t,u
te,t,e
to,t,o
na,n,a
ni,n,i
nu,n,u
ne,n,e
no,n,o
ha,h,a
hi,h,i
hu,h,u
he,h,e
ho,h,o
fa,f,a
fi,f,i
fu,f,u
fe,f,e
fo,f,o
ma,m,a
mi,m,i
mu,m,u
me,m,e
mo,m,o
ya,y,a
yi,y,i
yu,y,u
ye,y,e
yo,y,o
ra,r,a
ri,r,i
ru,r,u
re,r,e
ro,r,o
wa,w,a
wi,w,i
wu,w,u
we,w,e
wo,w,o
ga,g,a
gi,g,i
gu,g,u
ge,g,e
go,g,o
ja,j,a
ji,j,i
ju,j,u
je,j,e
jo,j,o
jya,j,a
jyi,j,i
jyu,j,u
jye,j,e
jyo,j,o
za,z,a
zi,z,i
zu,z,u
ze,z,e
zo,z,o
da,d,a
di,d,i
du,d,u
de,d,e
do,d,o
ba,b,a
bi,b,i
bu,b,u
be,b,e
bo,b,o
va,b,a
vi,b,i
vu,b,u
ve,b,e
vo,b,o
pa,p,a
pi,p,i
pu,p,u
pe,p,e
po,p,o
kya,ky,a
kyi,ky,i
kyu,ky,u
kye,ky,e
kyo,ky,o
tya,ty,a
tyi,ty,i
tyu,ty,u
tye,ty,e
tyo,ty,o
nya,ny,a
nyi,ny,i
nyu,ny,u
nye,ny,e
nyo,ny,o
hya,hy,a
hyi,hy,i
hyu,hy,u
hye,hy,e
hyo,hy,o
mya,my,a
myi,my,i
myu,my,u
mye,my,e
myo,my,o
rya,ry,a
ryi,ry,i
ryu,ry,u
rye,ry,e
ryo,ry,o
gya,gy,a
gyi,gy,i
gyu,gy,u
gye,gy,e
gyo,gy,o
bya,by,a
byi,by,i
byu,by,u
bye,by,e
byo,by,o
pya,py,a
pyi,py,i
pyu,py,u
pye,py,e
pyo,py,o
]]
romaji_to_hiragana = {}
hiragana_to_romaji = {}
phoneme_to_romaji = {}
config = {}
ext.readyn = function()
  repeat
    local input = ext.readstr()
    if input == "yes" or input == "YES" or input == "y" or input == "Y" then
      return true
    elseif input == "no" or input == "NO" or input == "n" or input == "N" then
      return false
    end
    ext.printstr("Please enter \"yes\" or \"no\".\n")
  until false
end

ext.promptyn = function(prompt)
  ext.printstr(prompt .. " (yes/y/no/n): ")
  return ext.readyn()
end


ext.listopt = function(options)
  repeat
    for i = 1, #options do
      ext.printstr("(" .. i .. ") " .. options[i] .. "\n")
    end
    ext.printstr("Please select an option from above: ")
    local choice = tonumber(ext.readstr())
    if choice ~= nil and 1 <= choice and choice <= #options then
      return choice
    end
  until false
end


utf8.decompose = function(str)
  local codepoints = {}
  for cp in utf8.codes(str) do
    codepoints[#codepoints + 1] = cp
  end
  codepoints[#codepoints + 1] = #str
  local result = {}
  for i = 1, #codepoints - 1 do
    result[#result + 1] = string.sub(str, codepoints[i], codepoints[i + 1] - 1)
  end
  return result
end


table.slice = function(t, start, finish)
  local result = {}
  for i = start, finish, 1 do
    result[i - start + 1] = t[i]
  end
  return result
end


table.append = function(dest, src)
  for i = 1, #src do
    dest[#dest + 1] = src[i]
  end
end


local function parse_csv_map(str)
  local lines = {{}}
  local line_start = 1
  if string.sub(str, #str, #str) ~= "\n" then
    str = str .. "\n"
  end
  for i = 1, #str do
    local ch = string.sub(str, i, i)
    if ch == "," or ch == "\n" or i == #str then
      local line = string.sub(str, line_start, i - 1)
      local current_line = lines[#lines]
      if string.byte(line, #line) == 13 then
        line = string.sub(line, 1, #line - 1)
      end
      current_line[#current_line + 1] = line
      if ch == "\n" or i == #str then
        lines[#lines + 1] = {}
      end
      line_start = i + 1
    end
  end
  return lines
end
local parsed = nil
parsed = parse_csv_map(hiragana_map_str)
for i = 1, #parsed do
  if #parsed[i] == 2 then
    local hiragana = parsed[i][1]
    local romaji = parsed[i][2]
    hiragana_to_romaji[hiragana] = romaji
    romaji_to_hiragana[romaji] = hiragana
  end
end
romaji_to_hiragana.ja = "じゃ"
romaji_to_hiragana.jo = "じょ"
romaji_to_hiragana.ju = "じゅ"
romaji_to_hiragana.wo = "を"
parsed = parse_csv_map(phoneme_map_str)
for i = 1, #parsed do
  if 1 < #parsed[i] then
    local phoneme_list = table.slice(parsed[i], 2, #parsed[i])
    local phoneme_key = table.concat(phoneme_list, ",")
    local phoneme = parsed[i][1]
    phoneme_to_romaji[phoneme_key] = phoneme
  end
end


-- These blank lines intentionally preserve the function line metadata in the original stripped Lua chunk.
-- They may originally have held formatting or source comments lost during compilation,
-- So do not add or remove lines here.








local function make_oto_entry(filename, name, offset, consonant_end, cutoff, preutt, overlap, end_time)
return { filename = filename, name = name,
offset = offset, consonant = consonant_end - offset, cutoff = end_time - cutoff,
preutt = preutt - offset, overlap = overlap - offset, }
end


local function oto_entry_to_params(entry)
  return math.round(entry.offset * 1000)
    .. "," .. math.round(entry.consonant * 1000)
    .. "," .. math.round(entry.cutoff * 1000)
    .. "," .. math.round(entry.preutt * 1000)
    .. "," .. math.round(entry.overlap * 1000)
end
local function oto_entry_to_line(entry)
return entry.filename .. "=" .. entry.name
.. config.suffix .. "," .. oto_entry_to_params(entry)
end


local function entries_to_string(entries)
  local lines = {}
  for i = 1, #entries do
    lines[i] = oto_entry_to_line(entries[i])
  end
  return table.concat(lines, "\n")
end


local function rename_duplicates(entries)
  local counts = {}
  for i = 1, #entries do
    local name = entries[i].name
    if counts[name] == nil then
      counts[name] = -1
    end
    counts[name] = counts[name] + 1
    if 1 <= counts[name] then
      entries[i].name = entries[i].name .. counts[name]
    end
  end
end

local function merge_existing_oto(in_oto, entry_ids, entries)
  local result = {}
  for i = 1, #entries do
    local key = entries[i].filename .. "=" .. entries[i].name
    if entry_ids[key] ~= nil then
      in_oto[entry_ids[key]].rest = oto_entry_to_params(entries[i])
    end
  end
end


local function is_silence(phoneme)
  return phoneme == "pau" or phoneme == "sil" or phoneme == "sil1" or phoneme == "sil2" or phoneme == "bre"




end


local function bool_to_string(value)
  if value == true then
    return "true"
  else
    return "false"
  end
end


local function get_phoneme_display(phoneme)
  if config.format.language == "japanese" then
    if phoneme == "N" then
      return "n"
    end
    if phoneme == "bre" then
      return "息"
    end
    if is_silence(phoneme) then
      return "-"
    end
    return phoneme
  else
    if is_silence(phoneme) then
      return "-"
    end
    return phoneme
  end
end


local function convert_grapheme(grapheme, target)
  if config.format.language == "japanese" then
    if config.format.variant == "hiragana" then
      if target == "romaji" then
        return hiragana_to_romaji[grapheme]
      end
      return grapheme
    end
    if config.format.variant == "romaji" then
      if target == "hiragana" then
        return romaji_to_hiragana[grapheme]
      end
      return grapheme
    end
    return grapheme
  end
  return grapheme




end

local function generate_cvvc_entries(ptr, label_phonemes)
  local entries = {}
  local name = label_phonemes.name
  local utt_end = label_phonemes[#label_phonemes].tend
  label_phonemes[#label_phonemes].grapheme = "-"
  for i = 2, #label_phonemes do
    if not is_silence(label_phonemes[i - 1].phoneme)
      or not is_silence(label_phonemes[i].phoneme)
    then
      local consonant_end = nil
      local vowel_start = nil
      local cutoff_point = nil
      local preutt_point = nil
      local offset_point = nil
      local prev_end = label_phonemes[i - 1].tbegin
      local curr_start = label_phonemes[i].tbegin
      local curr_end = label_phonemes[i].tend
      consonant_end = math.min(curr_start + 0.1, curr_start + (curr_end - curr_start) * 0.5)
      vowel_start = curr_end
      cutoff_point = curr_start
      preutt_point = math.max(curr_start - 0.01, prev_end + (curr_start - prev_end) * 0.5)
      offset_point = math.max(preutt_point - 0.1, prev_end + (curr_start - prev_end) * 0.3)
      if ptr.format.language == "english" or ext.is_consonant(ptr.ptr, label_phonemes[i].phoneme) then
        consonant_end = math.min(curr_start + 0.01, curr_start + (curr_end - curr_start) * 0.3)
      end
      local aliases = {}
      local phoneme_disp_prev = get_phoneme_display(label_phonemes[i - 1].phoneme)
      local phoneme_disp_curr = get_phoneme_display(label_phonemes[i].phoneme)
      local exclude_phoneme = false
      if config.format.language == "japanese" then
        local combined_key = label_phonemes[i - 1].phoneme .. "," .. label_phonemes[i].phoneme
        local combined_romaji = phoneme_to_romaji[combined_key]
        if combined_romaji ~= nil then
          if config.jpcvvc_romaji then
            aliases[#aliases + 1] = combined_romaji
            if not exclude_phoneme then
              exclude_phoneme = config.jpcvvc_exclusion
            end
          end
          if config.jpcvvc_hiragana and romaji_to_hiragana[combined_romaji] ~= nil then
            aliases[#aliases + 1] = romaji_to_hiragana[combined_romaji]
            if not exclude_phoneme then
              exclude_phoneme = config.jpcvvc_exclusion
            end
          end
        end
      end
      if config.cvvc_phoneme and not exclude_phoneme then
        aliases[#aliases + 1] = phoneme_disp_prev .. " " .. phoneme_disp_curr
      end
      for j = 1, #aliases do
        entries[#entries + 1] = make_oto_entry(name, aliases[j], offset_point, consonant_end, vowel_start, cutoff_point, preutt_point, utt_end)
      end
      if (is_silence(label_phonemes[i - 1].phoneme)
          or not ext.is_consonant(ptr.ptr, label_phonemes[i - 1].phoneme))
        and not ext.is_consonant(ptr.ptr, label_phonemes[i].phoneme)
        and not is_silence(label_phonemes[i].phoneme)
      then
        offset_point = cutoff_point
        cutoff_point = consonant_end - 0.01
        preutt_point = (cutoff_point + offset_point) / 2
        aliases = {}
        exclude_phoneme = false
        if config.format.language == "japanese" then
          if config.jpcvvc_hiragana and romaji_to_hiragana[phoneme_disp_curr] ~= nil then
            aliases[#aliases + 1] = romaji_to_hiragana[phoneme_disp_curr]
            if not exclude_phoneme then
              exclude_phoneme = config.jpcvvc_exclusion
            end
          end
        end
        if config.cvvc_phoneme and not exclude_phoneme then
          aliases[#aliases + 1] = phoneme_disp_curr
        end
        for j = 1, #aliases do
          entries[#entries + 1] = make_oto_entry(name, aliases[j], offset_point, consonant_end, vowel_start, cutoff_point, preutt_point, utt_end)
        end
      end
    end
  end
  return entries












end


local function generate_vcv_entries(ptr, label_phonemes)
  local entries = {}
  local name = label_phonemes.name
  local utt_end = label_phonemes[#label_phonemes].tend
  label_phonemes[#label_phonemes].grapheme = "-"
  for i = 2, #label_phonemes do
    if label_phonemes[i].grapheme ~= nil
      and (not is_silence(label_phonemes[i - 1].phoneme)
        or not is_silence(label_phonemes[i].phoneme))
    then
      local consonant_end = nil
      local vowel_start = nil
      local cutoff_point = nil
      local preutt_point = nil
      local offset_point = nil
      local prev_end = label_phonemes[i - 1].tbegin
      local curr_start = label_phonemes[i].tbegin
      local curr_end = label_phonemes[i].tend
      if ext.is_consonant(ptr.ptr, label_phonemes[i].phoneme) then
        local next_end = label_phonemes[i + 1].tend
        consonant_end = math.min(curr_end + 0.1, curr_end + (next_end - curr_end) * 0.2)
        vowel_start = math.max(next_end - 0.1, curr_end + (next_end - curr_end) * 0.66)
        cutoff_point = curr_end
        preutt_point = curr_start
        offset_point = math.max(0, preutt_point - 0.14)
      else
        consonant_end = math.min(curr_start + 0.1, curr_start + (curr_end - curr_start) * 0.2)
        vowel_start = math.max(curr_end - 0.1, curr_start + (curr_end - curr_start) * 0.66)
        cutoff_point = curr_start
        preutt_point = math.max(0, cutoff_point - 0.08)
        offset_point = math.max(0, cutoff_point - 0.16)
      end
      local aliases = {}
      local phoneme_disp_prev = get_phoneme_display(label_phonemes[i - 1].phoneme)
      local phoneme_disp_curr = get_phoneme_display(label_phonemes[i].phoneme)
      if config.vcv_phoneme then
        if ext.is_consonant(ptr.ptr, label_phonemes[i].phoneme) then
          local phoneme_disp_next = get_phoneme_display(label_phonemes[i + 1].phoneme)
          aliases[#aliases + 1] = phoneme_disp_prev .. " " .. phoneme_disp_curr .. " " .. phoneme_disp_next
        else
          aliases[#aliases + 1] = phoneme_disp_prev .. " " .. phoneme_disp_curr
        end
      end
      if config.format.language == "japanese" then
        local grapheme_romaji = convert_grapheme(label_phonemes[i].grapheme, "romaji")
        local grapheme_hiragana = convert_grapheme(label_phonemes[i].grapheme, "hiragana")
        if config.jpvcv_romaji and grapheme_romaji ~= nil then
          aliases[#aliases + 1] = phoneme_disp_prev .. " " .. grapheme_romaji
        end
        if config.jpvcv_hiragana and grapheme_hiragana ~= nil then
          aliases[#aliases + 1] = phoneme_disp_prev .. " " .. grapheme_hiragana
          if grapheme_hiragana == "お" then
            aliases[#aliases + 1] = phoneme_disp_prev .. " を"
          end
        end
      else
        aliases[#aliases + 1] = phoneme_disp_prev .. " " .. label_phonemes[i].grapheme
      end
      for j = 1, #aliases do
        entries[#entries + 1] = make_oto_entry(name, aliases[j], offset_point, consonant_end, vowel_start, cutoff_point, preutt_point, utt_end)
      end
    end
  end
  return entries






end















lmain = function(ptr)
  ext.printstr("Generating oto for language = "
    .. config.format.language
    .. ", variant = "
    .. config.format.variant
    .. "\n")
  local entries = {}
  for i = 1, #ptr.label_list do
    if config.generate_cvvc then
      table.append(entries, generate_cvvc_entries(ptr, ptr.label_list[i]))
    end
    if config.generate_vcv then
      table.append(entries, generate_vcv_entries(ptr, ptr.label_list[i]))
    end
  end
  if config.rename_dup then
    rename_duplicates(entries)
  end
  if config.entryids ~= nil then
    merge_existing_oto(config.in_oto, config.entryids, entries)
    local lines = {}
    for i = 1, #config.in_oto do
      lines[#lines + 1] = config.in_oto[i].first .. "," .. config.in_oto[i].rest
    end
    return (table.concat(lines, "\n"))
  else
    return (entries_to_string(entries))
  end



end

local function configure_japanese_vcv()
  ext.printstr("Available aliasing options for Japanese VCV voicebanks:\n")
  local choice = ext.listopt({
    "Hiragana only",
    "Romaji only",
    "Phonemes only",
    "Phonemes + Romaji",
    "Hiragana + Romaji",
    "Hiragana + Phonemes",
    "Hiragana + Romaji + Phonemes",
  })
  if choice == 1 or choice == 5 or choice == 6 or choice == 7 then
    config.jpvcv_hiragana = true
  end
  if choice == 2 or choice == 4 or choice == 5 or choice == 7 then
    config.jpvcv_romaji = true
  end
  if choice == 3 or choice == 4 or choice == 6 or choice == 7 then
    config.vcv_phoneme = true
  end
end

local function configure_japanese_cvvc()
ext.printstr("Available aliasing options for Japanese CVVC voicebanks:\n")
local choice = ext.listopt({
"Phonemes only",
"Phonemes + Romaji",
"Phonemes + Hiragana", "Phonemes + Hiragana + Romaji",
})
config.cvvc_phoneme = true
if choice == 2 or choice == 4 then
config.jpcvvc_romaji = true
end
if choice == 3 or choice == 4 then
config.jpcvvc_hiragana = true
end
if config.jpcvvc_romaji or config.jpcvvc_hiragana then config.jpcvvc_exclusion = ext.promptyn(
"Exclude phoneme alias for "
.. "those tokens whose romaji/hiragana alternatives exist?"
)
end
end

lconfig = function()
  config.generate_vcv = true
  config.generate_cvvc = true
  config.vcv_phoneme = false
  config.cvvc_phoneme = false
  config.jpvcv_romaji = false
  config.jpvcv_hiragana = false
  config.jpcvvc_romaji = false
  config.jpcvvc_hiragana = false
  config.jpcvvc_exclusion = false
  config.rename_dup = true
  config.suffix = ""
  ext.printstr("Configuring oto generation...\n")
  if config.in_oto ~= nil then
    local overwrite = ext.promptyn("oto.ini detected. Overwrite (yes) or update (no)?")
    if overwrite == false then
      local parsed = parse_csv_map(config.in_oto)
      config.entryids = {}
      config.in_oto = {}
      for i = 1, #parsed do
        local line = parsed[i][1]
        if line == nil then
          break
        end
        if #parsed[i] ~= 6 or string.find(line, "=") == nil then
          config.entryids = nil
          config.in_oto = nil
          ext.printstr("Warning: oto.ini is invalid. Switching to overwrite mode.\n")
          break
        end
        config.entryids[line] = i
        config.in_oto[i] = {
          first = line,
          rest = table.concat(table.slice(parsed[i], 2, 6), ","),
        }
      end
    end
  end
  if config.format.language == "japanese" then
    if config.entryids ~= nil then
      config.vcv_phoneme = true
      config.cvvc_phoneme = true
      config.jpvcv_romaji = true
      config.jpvcv_hiragana = true
      config.jpcvvc_romaji = true
      config.jpcvvc_hiragana = true
    else
      config.generate_vcv = ext.promptyn("Generate VCV entries?")
      config.generate_cvvc = ext.promptyn("Generate CVVC entries?")
      if config.generate_vcv then
        configure_japanese_vcv()
      end
      if config.generate_cvvc then
        configure_japanese_cvvc()
      end
    end
  elseif config.format.language == "english" then
    config.generate_vcv = false
    config.generate_cvvc = true
    config.cvvc_phoneme = true
  end
  if config.entryids ~= nil then
    config.rename_dup = false
  else
    config.rename_dup = ext.promptyn("Rename entries with duplicated name?")
  end
end
