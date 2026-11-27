"""Keep cumulative copies consistent with the corrected canonical sources."""
from pathlib import Path
import shutil
root=Path(__file__).resolve().parents[1]
common=['epoch.hpp','math.hpp','gravity.hpp','shgrav.hpp','mercury.hpp','mercury_parameters.hpp','eih.hpp','spice.hpp','rk8_table.hpp','integrator.hpp']
stages={1:common,2:common+['spline.hpp','stations.hpp','lighttime.hpp','trajectory.hpp'],3:common+['spline.hpp','stations.hpp','lighttime.hpp','trajectory.hpp','doppler.hpp'],4:common+['spline.hpp','stations.hpp','lighttime.hpp','trajectory.hpp','doppler.hpp','fit.hpp']}
for stage,directory in enumerate(sorted((root/'versions').iterdir()),1):
    for name in stages[stage]+['main.cpp']:
        shutil.copy2(root/'src'/name,directory/'src'/name)
    tests=directory/'tests';tests.mkdir(exist_ok=True)
    shutil.copy2(root/'tests/core_tests.cpp',tests/'core_tests.cpp')
    if stage>=2: shutil.copy2(root/'tests/physics_tests.cpp',tests/'physics_tests.cpp')
    if stage>=3: shutil.copy2(root/'tests/stage3_tests.cpp',tests/'stage3_tests.cpp')
    text=(root/'CMakeLists.txt').read_text(encoding='utf-8').replace('set(MESSENGER_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")','set(MESSENGER_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../..")').replace('set(MESSENGER_STAGE 4)',f'set(MESSENGER_STAGE {stage})')
    config=directory/'CMakeLists.txt'
    if not config.exists() or config.read_text(encoding='utf-8')!=text: config.write_text(text,encoding='utf-8')
print('Synchronized four cumulative versions')
