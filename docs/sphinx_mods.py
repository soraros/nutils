# Copyright (c) 2014 Evalf
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import inspect, pathlib, shutil, os, runpy, urllib.parse, shlex
import docutils.nodes, docutils.parsers.rst, docutils.statemachine
import sphinx.util.logging, sphinx.util.docutils, sphinx.addnodes
import nutils.log, nutils.matrix

project_root = pathlib.Path(__file__).parent.parent.resolve()

def process_signature(self, objtype, fullname, object, options, args, retann):
  if objtype in ('function', 'class', 'method'):
    signature = inspect.signature(object)
  else:
    return
  # Drop annotations from signature.
  signature = signature.replace(parameters=(param.replace(annotation=param.empty) for param in signature.parameters.values()),
                                return_annotation=inspect.Signature.empty)
  # Return a string representation of args and of the return annotation.  Note
  # that `str(signature)` would have included the return annotation if we
  # hadn't removed it above.
  return str(signature).replace('\\', '\\\\'), ''

def print_rst_autogen_header(*, file, src=None):
  print('..', file=file)
  print('   Automatically generated.  Edits are futile.', file=file)
  print(file=file)
  print(':autogenerated:', file=file)
  if src is not None:
    abssrc = src.resolve().relative_to(project_root)
    print(':autogeneratedfrom: {}'.format(abssrc), file=file)
  print(file=file)

def print_rst_h1(text, *, file):
  assert '\n' not in text
  print(file=file)
  print(text, file=file)
  print('='*len(text), file=file)
  print(file=file)

def print_rst_toctree(entries, *, relative_to, file, maxdepth=None):
  if relative_to.is_file():
    relative_to = relative_to.parent
  print(file=file)
  print('.. toctree::', file=file)
  if maxdepth is not None:
    print('   :maxdepth: {}'.format(maxdepth), file=file)
  print(file=file)
  for entry in entries:
    print('   {}'.format(entry.relative_to(relative_to).as_posix()), file=file)
  print(file=file)

def print_rst_label(name, *, file):
  print(file=file)
  print('.. _{}:'.format(name), file=file)
  print(file=file)

def copy_utime(src, dst):
  stat = os.stat(str(src))
  os.utime(str(dst), ns=(stat.st_atime_ns, stat.st_mtime_ns))

def generate_examples(app):
  logger = sphinx.util.logging.getLogger(__name__)

  dst_examples = pathlib.Path(app.srcdir)/'generated'/'examples'
  dst_examples.mkdir(parents=True, exist_ok=True)

  toc = []
  for src in sorted(project_root.glob('examples/*.py'), key=lambda item: str(item.with_suffix(''))):
    if src.name == '__init__.py':
      continue
    logger.info('generating examples... {}'.format(src))
    name = src.name
    dst = dst_examples/(name+'.rst')
    toc.append(dst)

    with src.open('r') as f_src, dst.open('w') as f_dst:
      print_rst_autogen_header(file=f_dst, src=src)
      # Add a label such that you can reference an example by
      # :ref:`examples/laplace.py`.
      print_rst_label('examples/{}'.format(name), file=f_dst)
      print_rst_h1(name, file=f_dst)
      print('.. exampledoc:: {}'.format(src.relative_to(project_root).as_posix()), file=f_dst)
    copy_utime(src, dst)

  logger.info('generating examples index...')
  index = dst_examples.with_suffix('.rst')
  with index.open('w') as f:
    print_rst_autogen_header(file=f)
    print_rst_h1('Examples', file=f)
    print_rst_toctree(toc, relative_to=index, maxdepth=1, file=f)

class LineIter:

  def __init__(self, lines):
    self._lines = iter(lines)
    self._index = -1
    self._next = None
    self.__next__()

  def __bool__(self):
    return self._next != StopIteration

  def __iter__(self):
    return self

  def __next__(self):
    if self._next == StopIteration:
      raise StopIteration
    value = self._index, self._next
    try:
      self._next = next(self._lines)
      self._index += 1
    except StopIteration:
      self._next = StopIteration
    return value

  @property
  def peek(self):
    if self._next == StopIteration:
      raise ValueError
    else:
      return self._next

class ExampleDocDirective(docutils.parsers.rst.Directive):

  has_content = False
  required_arguments = 1
  options_arguments = 0

  @staticmethod
  def _isdocline(line):
    line = line.lstrip()
    return line.rstrip() == '#' or line.startswith('# ')

  def run(self):
    logger = sphinx.util.logging.getLogger(__name__)
    nodes = []

    src = project_root/self.arguments[0]
    with src.open('r') as f:
      prevtype = None
      lines = LineIter(f)
      if lines and lines.peek.startswith('#!'):
        next(lines)
      while lines:
        if lines.peek.rstrip('\n') == '':
          next(lines)
        elif self._isdocline(lines.peek):
          # Collect all doc lines.
          contents = docutils.statemachine.ViewList()
          while lines and self._isdocline(lines.peek):
            i, line = next(lines)
            contents.append(line.lstrip()[2:], self.arguments[0], i)
          # Parse as rst into `node`.
          with sphinx.util.docutils.switch_source_input(self.state, contents):
            node = docutils.nodes.container()
            self.state.nested_parse(contents, 0, node)
          # Process sh roles.  Add links to logs.
          for sh_node in node.traverse(docutils.nodes.literal):
            if 'nutils_sh' not in sh_node:
              continue
            cmdline = sh_node.get('nutils_sh')
            cmdline_parts = tuple(shlex.split(cmdline))
            if cmdline_parts[:2] != ('python3', src.name):
              logger.warn('Not creating a log for {}.'.format(cmdline))
              continue
            log_link = sphinx.addnodes.only(expr='html')
            log_link.append(docutils.nodes.inline('', ' '))
            xref = sphinx.addnodes.pending_xref('', reftype='nutils-log', refdomain='std', reftarget=cmdline_parts[2:], script=src)
            xref += docutils.nodes.inline('', '(view log)', classes=['nutils-log-link'])
            log_link += xref
            sh_node.parent.insert(sh_node.parent.index(sh_node)+1, log_link)
          nodes.extend(node.children)
        else:
          # Collect all source lines.
          istart, line = next(lines)
          contents = [line]
          while lines and not self._isdocline(lines.peek):
            i, line = next(lines)
            contents.append(line)
          # Remove trailing empty lines.
          while contents and contents[-1].rstrip('\n') == '':
            del contents[-1]
          contents = ''.join(contents)
          # Create literal block.
          literal = docutils.nodes.literal_block(contents, contents)
          literal['language'] = 'python3'
          literal['linenos'] = True
          literal['highlight_args'] = dict(linenostart=istart+1)
          sphinx.util.nodes.set_source_info(self, literal)
          nodes.append(literal)

    return nodes

def role_sh(name, rawtext, text, lineno, inliner, options={}, context=[]):
  return [docutils.nodes.literal('', text, nutils_sh=text)], []

def create_log(app, env, node, contnode):
  logger = sphinx.util.logging.getLogger(__name__)

  if node['reftype'] == 'nutils-log':
    script = node.get('script')
    scriptname = str(script.relative_to(project_root))

    cmdline_args = node['reftarget']
    cmdline = ' '.join(map(shlex.quote, [scriptname, *cmdline_args]))

    target = '_logs/{}/index'.format(urllib.parse.quote(cmdline, safe='').replace('%', '+'))

    dst_log = (pathlib.Path(app.builder.outdir)/target).parent
    if dst_log.exists() and dst_log.stat().st_mtime > script.stat().st_mtime:
      logger.debug('Skip building log of {cmdline} because it already exists and '
                   'is newer than {script}.  Please touch {script} to force a rebuild.'
                   .format(script=scriptname, cmdline=cmdline))
    else:
      if dst_log.exists():
        logger.debug('purging old log files... {}'.format(dst_log))
        shutil.rmtree(str(dst_log))
      else:
        dst_log.parent.mkdir(parents=True, exist_ok=True)
      logger.info('creating log... {}'.format(cmdline))
      script_dict = runpy.run_path(str(script), run_name='__log_builder__')
      # Parse cmdline.
      params = inspect.signature(script_dict['main']).parameters.values()
      kwargs = {param.name: param.default for param in params}
      for arg in cmdline_args:
        if not arg:
          continue
        name, sep, value = arg.lstrip('-').partition('=')
        if not sep:
          value = not name.startswith('no')
          if not value:
            name = name[2:]
        if name not in kwargs:
          logger.error('unkown argument {!r}'.format(name))
          return
        default = kwargs[name]
        try:
          if isinstance(default, bool) and not isinstance(value, bool):
            raise Exception('boolean value should be specifiec as --{0}/--no{0}'.format(name))
          kwargs[name] = default.__class__(value)
        except Exception as e:
          logger.error('invalid argument for {!r}: {}'.format(name, e))
          return
      # Run script.
      with nutils.log.HtmlLog(str(dst_log)), nutils.matrix.backend('scipy'):
        script_dict['main'](**kwargs)
      (dst_log/'log.html').rename(dst_log/'index.html')

    refnode = docutils.nodes.reference('', '', internal=False, refuri=app.builder.get_relative_uri(env.docname, target))
    refnode.append(contnode)
    return refnode

def generate_api(app):
  logger = sphinx.util.logging.getLogger(__name__)

  nutils = project_root/'nutils'
  dst_root = pathlib.Path(app.srcdir)/'generated'/'api'
  dst_root.mkdir(parents=True, exist_ok=True)

  toc = []
  for src in sorted(nutils.glob('**/*.py')):
    if src == nutils/'__init__.py':
      continue
    module = '.'.join((src.parent if src.name == '__init__.py' else src.with_suffix('')).relative_to(nutils).parts)
    logger.info('generating api... {}'.format(module))
    dst = dst_root/(module+'.rst')
    toc.append(dst)
    with dst.open('w') as f:
      print_rst_autogen_header(file=f, src=src)
      print_rst_h1(module, file=f)
      print('.. automodule:: {}'.format('nutils.{}'.format(module)), file=f)
    copy_utime(src, dst)

  logger.info('generating api index...')
  index = dst_root.with_suffix('.rst')
  with index.open('w') as f:
    print_rst_autogen_header(file=f)
    print_rst_h1('API reference', file=f)
    print_rst_toctree(toc, relative_to=index, maxdepth=1, file=f)

def remove_generated(app, exception):
  logger = sphinx.util.logging.getLogger(__name__)
  generated = pathlib.Path(app.srcdir)/'generated'
  shutil.rmtree(str(generated), onerror=lambda f, p, e: logger.warning('failed to remove {}'.format(p)))

def setup(app):
  app.connect('autodoc-process-signature', process_signature)

  app.connect('builder-inited', generate_api)

  app.connect('builder-inited', generate_examples)
  app.add_directive('exampledoc', ExampleDocDirective)
  app.add_role('sh', role_sh)
  app.connect('missing-reference', create_log)

  app.connect('build-finished', remove_generated)

# vim: sts=2:sw=2:et
